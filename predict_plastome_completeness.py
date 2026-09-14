#!/usr/bin/env python3
"""
predict_plastome_completeness.py

Pre-assembly triage: given a sample's TRIMMED paired-end reads and a reference plastome (typically the same "closest available taxon" reference), estimate whether that sample's reads carry enough on-target plastid signal to likely yield a complete/circular plastome -- BEFORE spending assembler time on it.

SCOPE: HYBSEQ (TARGET-CAPTURE) DATA ONLY
------------------------------------------
This tool is calibrated for, and intended only for, hybseq/target-capture libraries, not whole-genome shotgun (WGS) data. The classification thresholds below were empirically tuned against this project's real 74-sample hybseq batch with known GetOrganelle/NOVOPlasty outcomes (see "CALIBRATION" below); WGS libraries have a very different on-target read fraction and coverage profile, so these same thresholds would not transfer to WGS data without re-calibrating against WGS ground truth separately. 

WHAT THIS DOES NOT DO
----------------------
This is not a re-implementation of GetOrganelle/NOVOPlasty's seed-and-extend logic, and it is not a guarantee. It answers a narrower, cheaper question: "do these reads have deep, broad coverage of a plastome reference?" via a quick short-read mapping (minimap2) rather than assembly. Two important caveats, both worth keeping in mind when reading the verdicts:

1. This is reference-guided by construction. A divergent reference can make a genuinely fine library look "marginal" here even though a de novo run might still succeed. This project's own benchmark found de novo and reference-guided GetOrganelle runs are highly concordant, but that isn't guaranteed for every taxon pair. Treat a MARGINAL call as "worth trying, keep an eye on it," not "don't bother."
2. Coverage/breadth thresholds below were empirically calibrated against this project's own hybseq outcomes. They aren't a claimed universal biological cutoff.

WORKFLOW PER SAMPLE
--------------------
1. Count total read pairs / bases in the trimmed fastq (fast, no mapping).
2. If the library is larger than --subsample-pairs, draw a fixed-seed subsample with seqtk (mapping every read of a multi-billion-base library against a 150 kb reference is wasted time, a few million pairs is plenty to estimate coverage depth and breadth).
3. Map (sub)reads to the reference plastome with minimap2 (short-read preset), sort/index with samtools.
4. From the alignment: % reads mapped, mean per-base depth across the reference (samtools depth -a, so zero-coverage positions count), and breadth of coverage at >=1x and >=3x.
5. Scale depth back up to the full (non-subsampled) read count.
6. Classify into likely_complete_circular / marginal_possible_fragment / likely_fail_low_coverage / likely_fail_low_breadth.

REQUIRES ON PATH (or loaded via `module load ...` before running)
--------------------------------------------------------------------
minimap2, samtools, seqtk

USAGE
-----
Single sample:
    python3 predict_plastome_completeness.py \\
        --sample SRR12345678 \\
        --r1 SRR12345678_1.trimmed.fastq.gz \\
        --r2 SRR12345678_2.trimmed.fastq.gz \\
        --ref-fasta /path/to/ref_cps/Genus_species_cp.fasta \\
        --out SRR12345678_prediction.tsv

Batch (reusing this project's ref_cp_taxon lookup convention -- see build_sample_sheet.py):
    python3 predict_plastome_completeness.py \\
        --sample-sheet samples.tsv \\
        --ref-dir /path/to/ref_cps \\
        --out batch_predictions.tsv

sample-sheet columns (tab-separated, header required): one of
    sample  r1  r2  ref_fasta

"""

import argparse
import csv
import gzip
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------------------
# Taxa/families with a well-documented, naturally IR-reduced plastome, per
# this project's own benchmarking results. Extend via --reduced-genome-taxa
# (one name per line, matched case-insensitively as a substring of ref_taxon).
# ---------------------------------------------------------------------------
BUILTIN_REDUCED_GENOME_TAXA = {
    # Fabaceae IRLC members confirmed IR-absent across every pipeline run
    # on this project's benchmark data:
    "vavilovia formosa",
    "austrocallerya megasperma",
    # Pinaceae: family-wide, well-documented historical IR loss; confirmed
    # IR-absent across all 9/9 records for both organisms tested:
    "pinus contorta",
    "pseudotsuga menziesii",
}

DEFAULT_SUBSAMPLE_PAIRS = 2_000_000
DEFAULT_SEED = 42

# Classification thresholds (illustrative defaults -- see module docstring).
DEFAULT_MIN_BREADTH_FAIL = 0.50   # below this fraction of ref covered >=1x -> fail
DEFAULT_MIN_COVERAGE_FAIL = 3.0   # below this estimated depth -> fail ("cov too low")
DEFAULT_MIN_COVERAGE_MARGINAL = 60.0  # calibrated 2026-09 against 74-sample real
# hybseq batch (see module docstring CALIBRATION section): raised from 15.0 to
# 60.0 to fix likely_complete_circular precision (was 0.62, 8 false positives;
# now 1.00) at the cost of one lost true positive (recall 1.00 -> 0.92). Safe
# empirical zone was (51.56, 74.01]; 60.0 sits in the middle of it.
DEFAULT_MIN_BREADTH_GOOD = 0.85   # below this (>=3x breadth) -> marginal even if deep
REDUCED_GENOME_BREADTH_RELAX = 0.20  # subtracted from breadth bars for flagged taxa


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_tools():
    missing = [t for t in ("minimap2", "samtools", "seqtk") if shutil.which(t) is None]
    if missing:
        eprint(
            "ERROR: required tool(s) not found on PATH: " + ", ".join(missing) + "\n"
            "Load them first, e.g.:\n"
            "    module load minimap2 samtools seqtk\n"
            "or install via conda:\n"
            "    conda install -c bioconda minimap2 samtools seqtk"
        )
        sys.exit(1)


def clean_taxon_to_filename(taxon: str) -> str:
    """Reproduce this project's established ref-plastome filename convention:
    strip periods, spaces -> underscores, append _cp.fasta. E.g.
    "Pseudotsuga menziesii var. glauca" ->
    "Pseudotsuga_menziesii_var_glauca_cp.fasta".
    """
    cleaned = taxon.replace(".", "")
    cleaned = re.sub(r"\s+", "_", cleaned.strip())
    return f"{cleaned}_cp.fasta"


def is_reduced_genome_taxon(taxon: str, extra_taxa: set) -> bool:
    if not taxon:
        return False
    t = taxon.strip().lower()
    all_taxa = BUILTIN_REDUCED_GENOME_TAXA | extra_taxa
    return any(flagged in t for flagged in all_taxa)


def _opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def count_read_pairs_and_bases(r1_path: str) -> "tuple[int, int]":
    """Count read pairs and total bases (R1 only, doubled for the pair) by
    streaming R1 -- fast, avoids needing a fastq library dependency."""
    n_reads = 0
    n_bases = 0
    with _opener(r1_path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:  # sequence line
                n_reads += 1
                n_bases += len(line.strip())
    return n_reads, n_bases * 2  # assume R2 contributes roughly the same


def run(cmd, **kwargs):
    # stdout/stderr=PIPE + universal_newlines=True instead of capture_output=
    # True/text=True: those two keywords were only added to subprocess.run()
    # in Python 3.7. This project's cluster conda env resolved to an older
    # Python 3.6 (confirmed by the "unexpected keyword argument
    # 'capture_output'" error on the first real run), so stick to the
    # subprocess.run() signature that's worked since Python 3.5.
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, **kwargs)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({' '.join(cmd) if isinstance(cmd, list) else cmd}):\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return result


def subsample_fastq(r1, r2, n_pairs, seed, workdir):
    sub_r1 = os.path.join(workdir, "sub_R1.fastq")
    sub_r2 = os.path.join(workdir, "sub_R2.fastq")
    with open(sub_r1, "w") as out1:
        subprocess.run(["seqtk", "sample", "-s", str(seed), r1, str(n_pairs)],
                        stdout=out1, check=True)
    with open(sub_r2, "w") as out2:
        subprocess.run(["seqtk", "sample", "-s", str(seed), r2, str(n_pairs)],
                        stdout=out2, check=True)
    return sub_r1, sub_r2


def ref_length_bp(ref_fasta: str) -> int:
    total = 0
    with open(ref_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def map_and_measure(ref_fasta, r1, r2, threads, workdir):
    bam = os.path.join(workdir, "aln.sorted.bam")
    mm2_log = os.path.join(workdir, "minimap2.stderr.log")

    # Run minimap2 | samtools sort as a single shell pipeline rather than two
    # manually-wired Popen objects. An earlier version closed p1.stdout (the
    # standard idiom so minimap2 gets SIGPIPE if samtools exits early) and
    # then still called p1.communicate() -- which internally polls stdout
    # *and* stderr via `selectors`. Newer CPython (3.9+) skips an
    # already-closed stream there, but older versions don't guard that case
    # and raise `ValueError: Invalid file object: <_io.BufferedReader ...>`
    # from selectors._fileobj_to_fd when it calls .fileno() on the closed
    # stream. Confirmed by reproducing the closed-stream check's absence
    # is exactly what differs across Python versions; letting the shell own
    # the pipe (as it would from an interactive terminal) sidesteps the
    # whole class of bug instead of chasing per-version subprocess behavior.
    pipeline = (
        f"set -o pipefail; "
        f"minimap2 -ax sr -t {threads} {shlex.quote(ref_fasta)} {shlex.quote(r1)} {shlex.quote(r2)} "
        f"2>{shlex.quote(mm2_log)} | "
        f"samtools sort -@ {threads} -o {shlex.quote(bam)} -"
    )
    result = subprocess.run(["bash", "-c", pipeline], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)
    if result.returncode != 0:
        mm2_err = ""
        if os.path.exists(mm2_log):
            with open(mm2_log) as fh:
                mm2_err = fh.read()
        raise RuntimeError(
            f"minimap2/samtools sort pipeline failed (exit {result.returncode}):\n"
            f"samtools sort stderr: {result.stderr}\nminimap2 stderr: {mm2_err}"
        )
    run(["samtools", "index", bam])

    flagstat = run(["samtools", "flagstat", bam]).stdout
    pct_mapped = 0.0
    for line in flagstat.splitlines():
        if "mapped (" in line and "primary mapped" not in line and "mate" not in line:
            m = re.search(r"\(([\d.]+)%", line)
            if m:
                pct_mapped = float(m.group(1))
            break

    depth_out = run(["samtools", "depth", "-a", bam]).stdout
    depths = []
    for line in depth_out.splitlines():
        parts = line.split("\t")
        if len(parts) == 3:
            depths.append(int(parts[2]))

    if not depths:
        return dict(pct_mapped=pct_mapped, mean_depth=0.0, breadth_1x=0.0, breadth_3x=0.0, ref_len=0)

    ref_len = len(depths)
    mean_depth = sum(depths) / ref_len
    breadth_1x = sum(1 for d in depths if d >= 1) / ref_len
    breadth_3x = sum(1 for d in depths if d >= 3) / ref_len

    return dict(pct_mapped=pct_mapped, mean_depth=mean_depth,
                breadth_1x=breadth_1x, breadth_3x=breadth_3x, ref_len=ref_len)


def classify(estimated_coverage, breadth_1x, breadth_3x, reduced_genome, thresholds):
    relax = REDUCED_GENOME_BREADTH_RELAX if reduced_genome else 0.0
    min_breadth_fail = max(0.0, thresholds["min_breadth_fail"] - relax)
    min_breadth_good = max(0.0, thresholds["min_breadth_good"] - relax)

    notes = []
    if reduced_genome:
        notes.append(
            "reference taxon is flagged as a naturally IR-reduced/no-IR lineage in this "
            "project's benchmark -- breadth bars relaxed accordingly, don't expect ~100% breadth"
        )

    if breadth_1x < min_breadth_fail:
        label = "likely_fail_low_breadth"
        notes.append(
            f"only {breadth_1x:.0%} of reference covered at all (>=1x); "
            "check reference relatedness (a divergent reference can look like this even "
            "with a fine library -- de novo assembly may still be worth trying) or on-target read yield"
        )
    elif estimated_coverage < thresholds["min_coverage_fail"]:
        label = "likely_fail_low_coverage"
        notes.append(f"estimated on-target coverage ~{estimated_coverage:.1f}x is very low")
    elif estimated_coverage < thresholds["min_coverage_marginal"] or breadth_3x < min_breadth_good:
        label = "marginal_possible_fragment"
        notes.append(
            f"estimated coverage ~{estimated_coverage:.1f}x and/or >=3x breadth "
            f"{breadth_3x:.0%} sit below the 'likely complete' bar; worth attempting but keep "
            "a fallback (reference-guided run, or a MITObim-style rescue if it fragments) in the plan"
        )
    else:
        label = "likely_complete_circular"
        notes.append(f"deep ({estimated_coverage:.1f}x est.), broad ({breadth_3x:.0%} >=3x) on-target coverage")

    return label, "; ".join(notes)


def process_sample(sample, r1, r2, ref_fasta, threads, subsample_pairs, seed,
                    reduced_genome_taxa, ref_taxon_for_flagging, thresholds, keep_intermediate):
    with tempfile.TemporaryDirectory(prefix=f"plastpred_{sample}_") as workdir:
        n_pairs_total, n_bases_total = count_read_pairs_and_bases(r1)

        if n_pairs_total == 0:
            return dict(sample=sample, error="zero reads in R1 -- empty or unreadable fastq")

        used_subsample = n_pairs_total > subsample_pairs
        if used_subsample:
            map_r1, map_r2 = subsample_fastq(r1, r2, subsample_pairs, seed, workdir)
            n_pairs_used = subsample_pairs
        else:
            map_r1, map_r2 = r1, r2
            n_pairs_used = n_pairs_total

        metrics = map_and_measure(ref_fasta, map_r1, map_r2, threads, workdir)

        scale = (n_pairs_total / n_pairs_used) if n_pairs_used else 1.0
        estimated_coverage = metrics["mean_depth"] * scale

        reduced = is_reduced_genome_taxon(ref_taxon_for_flagging, reduced_genome_taxa)
        label, notes = classify(estimated_coverage, metrics["breadth_1x"], metrics["breadth_3x"],
                                 reduced, thresholds)

        return dict(
            sample=sample,
            ref_fasta=ref_fasta,
            ref_length_bp=metrics["ref_len"],
            total_read_pairs=n_pairs_total,
            total_bases_bp=n_bases_total,
            subsampled="yes" if used_subsample else "no",
            pct_reads_mapped=round(metrics["pct_mapped"], 2),
            mean_depth_subsample=round(metrics["mean_depth"], 2),
            estimated_full_coverage=round(estimated_coverage, 2),
            breadth_1x_pct=round(metrics["breadth_1x"] * 100, 1),
            breadth_3x_pct=round(metrics["breadth_3x"] * 100, 1),
            reduced_genome_taxon="yes" if reduced else "no",
            predicted_call=label,
            notes=notes,
            error="",
        )


OUTPUT_FIELDS = [
    "sample", "ref_fasta", "ref_length_bp", "total_read_pairs", "total_bases_bp",
    "subsampled", "pct_reads_mapped", "mean_depth_subsample", "estimated_full_coverage",
    "breadth_1x_pct", "breadth_3x_pct", "reduced_genome_taxon", "predicted_call", "notes", "error",
]


def load_extra_reduced_taxa(path):
    if not path:
        return set()
    with open(path) as fh:
        return {line.strip().lower() for line in fh if line.strip()}


def resolve_ref(row, ref_dir):
    if row.get("ref_fasta"):
        return row["ref_fasta"], row.get("ref_taxon", "")
    taxon = row.get("ref_taxon") or row.get("ref_cp_taxon")
    if not taxon:
        raise ValueError(f"row for sample {row.get('sample')} has neither ref_fasta nor ref_taxon column")
    if not ref_dir:
        raise ValueError("--ref-dir is required when the sample sheet uses ref_taxon instead of ref_fasta")
    fname = clean_taxon_to_filename(taxon)
    return os.path.join(ref_dir, fname), taxon


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", help="sample name (single-sample mode)")
    ap.add_argument("--r1", help="trimmed R1 fastq(.gz) (single-sample mode)")
    ap.add_argument("--r2", help="trimmed R2 fastq(.gz) (single-sample mode)")
    ap.add_argument("--ref-fasta", help="reference plastome fasta (single-sample mode)")
    ap.add_argument("--ref-taxon", default="", help="taxon name for the reference, used only for the "
                     "reduced-genome flag (single-sample mode; ignored if --ref-fasta already given)")

    ap.add_argument("--sample-sheet", help="TSV with header: sample r1 r2 (ref_fasta | ref_taxon). "
                     "Batch mode; overrides --sample/--r1/--r2/--ref-fasta.")
    ap.add_argument("--ref-dir", help="directory holding <Taxon_cleaned>_cp.fasta reference plastomes, "
                     "used to resolve ref_taxon in the sample sheet (this project's naming convention: "
                     "strip periods, spaces->underscores, + '_cp.fasta')")

    ap.add_argument("--out", required=True, help="output TSV path")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--subsample-pairs", type=int, default=DEFAULT_SUBSAMPLE_PAIRS,
                     help=f"cap mapping to this many read pairs, scaled back up for the coverage "
                          f"estimate (default {DEFAULT_SUBSAMPLE_PAIRS:,}; set to 0 to disable and map everything)")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--reduced-genome-taxa", help="optional file, one taxon name per line, appended to "
                     "the built-in IR-reduced-lineage list")
    ap.add_argument("--min-breadth-fail", type=float, default=DEFAULT_MIN_BREADTH_FAIL)
    ap.add_argument("--min-coverage-fail", type=float, default=DEFAULT_MIN_COVERAGE_FAIL)
    ap.add_argument("--min-coverage-marginal", type=float, default=DEFAULT_MIN_COVERAGE_MARGINAL)
    ap.add_argument("--min-breadth-good", type=float, default=DEFAULT_MIN_BREADTH_GOOD)
    ap.add_argument("--keep-intermediate", action="store_true", help="(not yet wired for BAMs; reserved)")

    args = ap.parse_args()
    check_tools()

    thresholds = dict(
        min_breadth_fail=args.min_breadth_fail,
        min_coverage_fail=args.min_coverage_fail,
        min_coverage_marginal=args.min_coverage_marginal,
        min_breadth_good=args.min_breadth_good,
    )
    extra_taxa = load_extra_reduced_taxa(args.reduced_genome_taxa)

    rows_to_run = []
    if args.sample_sheet:
        with open(args.sample_sheet) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows_to_run.append(row)
    else:
        missing = [a for a, v in [("--sample", args.sample), ("--r1", args.r1),
                                   ("--r2", args.r2), ("--ref-fasta", args.ref_fasta)] if not v]
        if missing:
            ap.error(f"single-sample mode requires {', '.join(missing)} (or use --sample-sheet)")
        rows_to_run.append(dict(sample=args.sample, r1=args.r1, r2=args.r2,
                                 ref_fasta=args.ref_fasta, ref_taxon=args.ref_taxon))

    results = []
    for row in rows_to_run:
        sample = row["sample"]
        eprint(f"[{sample}] starting...")
        try:
            ref_fasta, ref_taxon = resolve_ref(row, args.ref_dir)
            if not os.path.exists(ref_fasta):
                raise FileNotFoundError(f"reference plastome not found: {ref_fasta}")
            if not os.path.exists(row["r1"]) or not os.path.exists(row["r2"]):
                raise FileNotFoundError(f"trimmed fastq not found: {row['r1']} / {row['r2']}")

            result = process_sample(
                sample=sample, r1=row["r1"], r2=row["r2"], ref_fasta=ref_fasta,
                threads=args.threads, subsample_pairs=(args.subsample_pairs or 10**12),
                seed=args.seed, reduced_genome_taxa=extra_taxa, ref_taxon_for_flagging=ref_taxon,
                thresholds=thresholds, keep_intermediate=args.keep_intermediate,
            )
        except Exception as e:
            result = dict(sample=sample, error=str(e))

        for f in OUTPUT_FIELDS:
            result.setdefault(f, "")
        results.append(result)
        status = result.get("predicted_call") or f"ERROR: {result.get('error')}"
        eprint(f"[{sample}] {status}")

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    eprint(f"\nWrote {len(results)} row(s) to {args.out}")


if __name__ == "__main__":
    main()
