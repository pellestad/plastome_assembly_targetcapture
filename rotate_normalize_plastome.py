#!/usr/bin/env python3
"""
rotate_normalize_plastome.py

Set a uniform starting point and orientation for complete, circular plastid genome assemblies so that every sample in a dataset is directly comparable (same start coordinate, same strand, same SSC state).

Two steps, run per assembly:

  1. ROTATE to trnH-GUG. BLAST the trnH-GUG gene against the assembly, reverse-complement the whole molecule if the hit is on the minus strand (so the gene ends up read 5'->3'), then rotate the circle so the first base of the gene is position 1.

  2. NORMALIZE SSC ORIENTATION using ndhF. Self-BLAST the rotated molecule to find the inverted-repeat (IR) pair, use IR coordinates to split the molecule into LSC and SSC, then BLAST ndhF against the SSC region only. If ndhF is not on the expected ("canonical") strand, reverse-complement just the SSC segment in place (LSC/IR content and the trnH start point are untouched).

Requires NCBI BLAST+ (blastn) on PATH and Biopython.

USAGE
-----
Single sample:
    conda activate biopython   # or any env with biopython + blastn on PATH

    python rotate_normalize_plastome.py \
        --query-trnh /path/to/trnH_GUG_datasets/ncbi_dataset/data/gene.fna \
        --query-ndhf /path/to/ndhF_datasets/ncbi_dataset/data/gene.fna \
        --input sample1.complete.path_sequence.fasta \
        --outdir rotated_plastomes/

Batch over a directory (or list) of samples:
    python rotate_normalize_plastome.py \
        --query-trnh /path/to/trnH_GUG_datasets/ncbi_dataset/data/gene.fna \
        --query-ndhf /path/to/ndhF_datasets/ncbi_dataset/data/gene.fna \
        --input-dir /path/to/complete_plastomes/ \
        --pattern "*_complete_plastome.fasta" \
        --recursive \
        --outdir rotated_plastomes/ \
        2>&1 | tee rotate_normalize_run.log
"""

import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq

BLASTN_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand",
]


class SampleError(Exception):
    """Raised for a single-sample failure; caught in the batch loop so one
    bad sample does not stop the rest."""


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def write_fasta(path, header, seq, wrap=70):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


def revcomp(seq):
    return str(Seq(seq).reverse_complement())


def run_blastn(query_fasta, subject_fasta, evalue="1e-10", extra_args=None):
    # IMPORTANT: the `blastn` executable's default -task is actually
    # "megablast" (word_size 28), tuned for near-identical sequences. That
    # silently fails to seed any alignment once identity drops into the
    # ~75-90% range that's completely normal for a cross-species gene
    # comparison (e.g. ndhF against a divergent ortholog) -- it looks like
    # "no hit" rather than an error. Force the traditional blastn task
    # (word_size 11) everywhere so moderate divergence is actually found;
    # specificity is still controlled by the evalue/pident/length filters
    # applied by callers, not by the seeding word size.
    cmd = [
        "blastn",
        "-task", "blastn",
        "-query", query_fasta,
        "-subject", subject_fasta,
        "-outfmt", "6 " + " ".join(BLASTN_FIELDS),
        "-evalue", evalue,
    ]
    if extra_args:
        cmd += extra_args
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise SampleError(f"blastn failed: {result.stderr.strip()}")
    hits = []
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        hits.append(dict(zip(BLASTN_FIELDS, line.split("\t"))))
    return hits


def best_hit(hits, strand=None):
    """Highest-bitscore hit, optionally restricted to one strand. Filtering
    by strand BEFORE taking the max (rather than taking the unconditional
    max and checking its strand afterward) matters specifically in
    rotate_to_trnh's post-reverse-complement re-search: if two loci score
    at or near a tie against the query (e.g. a genuinely duplicated
    trnH-GUG copy sitting inside the IR itself, which happens in some
    lineages when the LSC/IR boundary falls just past the gene -- the two
    IR copies are near-identical, so both score similarly against the
    query), an unconditional max can end up pointing at whichever of the
    two tied loci happens to sort first, on either strand, rather than at
    the specific locus already confirmed as the anchor by the first
    search."""
    if strand is not None:
        hits = [h for h in hits if h["sstrand"] == strand]
    if not hits:
        return None
    return max(hits, key=lambda h: float(h["bitscore"]))


# --------------------------------------------------------------------------
# step 1: rotate to trnH-GUG
# --------------------------------------------------------------------------

def rotate_to_trnh(seq, query_trnh, tmp_dir, sample, min_pident, min_len):
    subj = os.path.join(tmp_dir, f"{sample}.trnh_subject.fasta")
    write_fasta(subj, sample, seq)
    hits = [h for h in run_blastn(query_trnh, subj)
            if float(h["pident"]) >= min_pident and int(h["length"]) >= min_len]
    hit = best_hit(hits)
    if hit is None:
        raise SampleError("no trnH-GUG hit found against this assembly")

    if hit["sstrand"] == "minus":
        # Reverse-complement the whole molecule so the gene reads 5'->3',
        # then re-BLAST against the flipped sequence to get clean plus-
        # strand coordinates rather than doing the coordinate arithmetic
        # by hand.
        seq = revcomp(seq)
        write_fasta(subj, sample, seq)
        hits2 = [h for h in run_blastn(query_trnh, subj)
                 if float(h["pident"]) >= min_pident and int(h["length"]) >= min_len]
        # Restrict to plus-strand hits before ranking, rather than taking
        # the unconditional best hit and complaining if it isn't plus --
        # see best_hit()'s docstring for why the latter is fragile to
        # near-tied competing loci (e.g. a real trnH-GUG duplicate inside
        # the IR).
        hit2 = best_hit(hits2, strand="plus")
        if hit2 is None:
            raise SampleError(
                "trnH-GUG hit lost after reverse-complementing (no plus-strand "
                "hit found even though a minus-strand hit existed before)"
            )
        anchor = int(hit2["sstart"])
        orig_strand = "minus"
    else:
        anchor = int(hit["sstart"])
        orig_strand = "plus"

    rotated = seq[anchor - 1:] + seq[:anchor - 1]
    return rotated, anchor, orig_strand


# --------------------------------------------------------------------------
# step 2: locate the IR pair, then normalize SSC orientation via ndhF
# --------------------------------------------------------------------------

def find_ir_pair(seq, tmp_dir, sample, min_ir_len, min_pident):
    subj = os.path.join(tmp_dir, f"{sample}.self_subject.fasta")
    write_fasta(subj, sample, seq)
    hits = run_blastn(subj, subj, evalue="1e-50")

    candidates = []
    for h in hits:
        if h["sstrand"] != "minus":
            continue  # IRs are inverted repeats of each other -> opposite strand
        length = int(h["length"])
        pident = float(h["pident"])
        if length < min_ir_len or pident < min_pident:
            continue
        qstart, qend = int(h["qstart"]), int(h["qend"])
        sstart, send = int(h["sstart"]), int(h["send"])
        copy1 = (min(qstart, qend), max(qstart, qend))
        copy2 = (min(sstart, send), max(sstart, send))
        if copy1[0] > copy2[0]:
            copy1, copy2 = copy2, copy1
        if copy1[1] >= copy2[0]:
            continue  # overlapping hit, not two distinct IR copies
        candidates.append((length, pident, copy1, copy2))

    if not candidates:
        return None

    candidates.sort(key=lambda c: c[0], reverse=True)
    length, pident, copy1, copy2 = candidates[0]
    return {"copy1": copy1, "copy2": copy2, "length": length, "pident": pident}


def diagnose_ir_wholegenome(seq, tmp_dir, sample):
    """Only called when no self-BLAST hit cleared the IR thresholds. Reruns
    the self-BLAST at a relaxed evalue (the primary search uses a strict
    1e-50 to keep noise out of real IR detection) purely to report the
    best *sub-threshold* opposite-strand, non-overlapping self-hit found,
    if any -- so a warning can say whether there's a weak/short IR just
    under the cutoff (worth relaxing --min-ir-len/--ir-min-pident) or
    genuinely nothing repeat-like in the molecule (consistent with real
    IR loss -- some plant lineages, e.g. the legume IR-lacking clade,
    genuinely lack one or both IR copies)."""
    subj = os.path.join(tmp_dir, f"{sample}.self_subject.fasta")
    hits = run_blastn(subj, subj, evalue="1e-5")
    best = None
    for h in hits:
        if h["sstrand"] != "minus":
            continue
        qstart, qend = int(h["qstart"]), int(h["qend"])
        sstart, send = int(h["sstart"]), int(h["send"])
        lo1, hi1 = min(qstart, qend), max(qstart, qend)
        lo2, hi2 = min(sstart, send), max(sstart, send)
        if lo1 > lo2:
            lo1, hi1, lo2, hi2 = lo2, hi2, lo1, hi1
        if hi1 >= lo2:
            continue  # overlapping, not two distinct copies
        length = int(h["length"])
        pident = float(h["pident"])
        if best is None or length > best["length"]:
            best = {"length": length, "pident": pident}
    if best is None:
        return {"best_len": "NA", "best_pident": "NA"}
    return {"best_len": best["length"], "best_pident": f"{best['pident']:.2f}"}


def normalize_ssc(seq, query_ndhf, ir_info, tmp_dir, sample,
                   canonical_strand, min_pident, min_len):
    copy1, copy2 = ir_info["copy1"], ir_info["copy2"]
    ssc_start, ssc_end = copy1[1] + 1, copy2[0] - 1  # 1-based, inclusive
    ssc_bounds = (ssc_start, ssc_end)
    if ssc_end <= ssc_start:
        return seq, None, False, "SSC region between IR copies is empty/invalid", ssc_bounds

    ssc_seq = seq[ssc_start - 1:ssc_end]
    subj = os.path.join(tmp_dir, f"{sample}.ssc_subject.fasta")
    write_fasta(subj, sample, ssc_seq)
    hits = [h for h in run_blastn(query_ndhf, subj)
            if float(h["pident"]) >= min_pident and int(h["length"]) >= min_len]
    hit = best_hit(hits)
    if hit is None:
        return seq, None, False, "ndhF not found within the identified SSC region", ssc_bounds

    strand = hit["sstrand"]
    if strand == canonical_strand:
        return seq, strand, False, None, ssc_bounds

    flipped_ssc = revcomp(ssc_seq)
    new_seq = seq[:ssc_start - 1] + flipped_ssc + seq[ssc_end:]
    return new_seq, strand, True, None, ssc_bounds


def diagnose_ndhf_wholegenome(seq, query_ndhf, ssc_bounds, tmp_dir, sample):
    """Only called when ndhF wasn't found inside the SSC slice. BLASTs ndhF
    against the *whole* rotated molecule with relaxed thresholds, purely to
    help distinguish 'gene present but outside the computed SSC window'
    (miscalculated SSC bounds) from 'gene genuinely absent/too divergent'
    (real biological ndh-suite loss, or a mismatched query file) -- it does
    not affect the normalization decision itself."""
    subj = os.path.join(tmp_dir, f"{sample}.wholegenome_subject.fasta")
    write_fasta(subj, sample, seq)
    # length floor of 100bp keeps this diagnostic meaningful -- with the
    # standard blastn task's small word size and a relaxed evalue, a much
    # shorter floor (e.g. 15-20bp) starts returning short, low-information
    # matches that occur by chance even in unrelated sequence.
    hits = [h for h in run_blastn(query_ndhf, subj, evalue="1")
            if int(h["length"]) >= 100]
    hit = best_hit(hits)
    if hit is None:
        return {"best_pident": "NA", "best_len": "NA", "best_start": "NA",
                "best_end": "NA", "best_strand": "NA", "in_ssc": "NA"}
    sstart, send = int(hit["sstart"]), int(hit["send"])
    lo, hi = min(sstart, send), max(sstart, send)
    in_ssc = ssc_bounds is not None and ssc_bounds[0] <= lo and hi <= ssc_bounds[1]
    return {
        "best_pident": f"{float(hit['pident']):.1f}",
        "best_len": hit["length"],
        "best_start": lo,
        "best_end": hi,
        "best_strand": hit["sstrand"],
        "in_ssc": in_ssc,
    }


# --------------------------------------------------------------------------
# per-sample driver
# --------------------------------------------------------------------------

def process_one(path, args, tmp_dir):
    sample = os.path.splitext(os.path.basename(path))[0]
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise SampleError(
            f"expected exactly 1 sequence in a complete plastome fasta, found {len(records)}"
        )
    seq = str(records[0].seq).upper()

    rotated, anchor, orig_strand = rotate_to_trnh(
        seq, args.query_trnh, tmp_dir, sample, args.min_pident, args.min_len
    )

    warnings = []
    ir_info = find_ir_pair(rotated, tmp_dir, sample, args.min_ir_len, args.ir_min_pident)

    ndhf_strand = None
    ssc_flipped = False
    final_seq = rotated
    ssc_bounds = None
    diag = None

    ir_diag = None
    if ir_info is None:
        ir_diag = diagnose_ir_wholegenome(rotated, tmp_dir, sample)
        if ir_diag["best_len"] == "NA":
            warnings.append(
                "no inverted-repeat pair found via self-BLAST, even at a relaxed "
                "evalue; SSC normalization skipped -- consistent with real IR loss "
                "in this lineage (e.g. the legume IR-lacking clade), or a "
                "genuinely IR-less/degraded assembly"
            )
        else:
            warnings.append(
                f"no inverted-repeat pair cleared the thresholds (--min-ir-len "
                f"{args.min_ir_len}, --ir-min-pident {args.ir_min_pident}), but the "
                f"best sub-threshold self-hit was {ir_diag['best_len']}bp at "
                f"{ir_diag['best_pident']}% identity; SSC normalization skipped -- "
                "if that's a real (short/divergent) IR for this lineage, rerun with "
                "relaxed --min-ir-len/--ir-min-pident"
            )
    else:
        L = len(rotated)
        len_between = ir_info["copy2"][0] - ir_info["copy1"][1] - 1
        len_wrap = (L - ir_info["copy2"][1]) + (ir_info["copy1"][0] - 1)
        if len_wrap <= len_between:
            warnings.append(
                f"unexpected IR geometry (between-IR region {len_between} bp >= "
                f"wrap-around region {len_wrap} bp) -- trnH-GUG may not sit in the "
                "expected LSC position for this sample; SSC normalization skipped, "
                "recommend manual review"
            )
        else:
            final_seq, ndhf_strand, ssc_flipped, note, ssc_bounds = normalize_ssc(
                rotated, args.query_ndhf, ir_info, tmp_dir, sample,
                args.canonical_ndhf_strand, args.ndhf_min_pident, args.ndhf_min_len,
            )
            if note:
                warnings.append(note)
                if "not found within the identified SSC region" in note:
                    diag = diagnose_ndhf_wholegenome(
                        rotated, args.query_ndhf, ssc_bounds, tmp_dir, sample
                    )
                    if diag["best_start"] == "NA":
                        warnings.append(
                            "diagnostic: ndhF not found anywhere in the whole molecule "
                            "even at relaxed thresholds -- check the ndhF query fasta, "
                            "or this taxon may genuinely lack/have a highly divergent ndhF"
                        )
                    else:
                        loc = "inside" if diag["in_ssc"] else "OUTSIDE"
                        warnings.append(
                            f"diagnostic: best whole-genome ndhF hit is {loc} the computed "
                            f"SSC window (pident={diag['best_pident']}%, len={diag['best_len']}, "
                            f"pos={diag['best_start']}-{diag['best_end']}, strand={diag['best_strand']})"
                        )

    out_path = os.path.join(args.outdir, f"{sample}.rotated.fasta")
    write_fasta(
        out_path,
        f"{sample} start=trnH-GUG SSC_normalized={'yes' if ssc_flipped or ndhf_strand == args.canonical_ndhf_strand else 'unverified'}",
        final_seq,
    )

    return {
        "sample": sample,
        "status": "OK",
        "length": len(final_seq),
        "trnH_anchor_pos_before_rotation": anchor,
        "trnH_orig_strand": orig_strand,
        "ir_len": ir_info["length"] if ir_info else "NA",
        "ir_pident": f"{ir_info['pident']:.2f}" if ir_info else "NA",
        "ir_diag_best_len": ir_diag["best_len"] if ir_diag else "NA",
        "ir_diag_best_pident": ir_diag["best_pident"] if ir_diag else "NA",
        "ssc_start": ssc_bounds[0] if ssc_bounds else "NA",
        "ssc_end": ssc_bounds[1] if ssc_bounds else "NA",
        "ssc_len": (ssc_bounds[1] - ssc_bounds[0] + 1) if ssc_bounds else "NA",
        "ndhF_strand_observed": ndhf_strand or "NA",
        "ssc_flipped": ssc_flipped,
        "ndhF_diag_best_pident": diag["best_pident"] if diag else "NA",
        "ndhF_diag_best_pos": f"{diag['best_start']}-{diag['best_end']}" if diag and diag["best_start"] != "NA" else "NA",
        "ndhF_diag_best_strand": diag["best_strand"] if diag else "NA",
        "ndhF_diag_in_ssc": diag["in_ssc"] if diag else "NA",
        "warnings": "; ".join(warnings),
        "output": out_path,
    }


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def collect_inputs(args):
    if args.input:
        return [args.input]
    pattern = os.path.join(args.input_dir, "**", args.pattern) if args.recursive \
        else os.path.join(args.input_dir, args.pattern)
    files = sorted(glob.glob(pattern, recursive=args.recursive))
    return files


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--query-trnh", required=True, help="trnH-GUG gene fasta (anchor for rotation)")
    p.add_argument("--query-ndhf", required=True, help="ndhF gene fasta (anchor for SSC orientation)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--input", help="single complete plastome fasta")
    src.add_argument("--input-dir", help="directory to search for plastome fastas")
    p.add_argument("--pattern", default="*.fasta", help="glob pattern used with --input-dir (default: *.fasta)")
    p.add_argument("--recursive", action="store_true", help="search --input-dir recursively")
    p.add_argument("--outdir", required=True, help="directory for rotated/normalized output fastas + log")
    p.add_argument("--canonical-ndhf-strand", choices=["plus", "minus"], default="minus",
                    help="strand ndhF should end up on after normalization (VERIFY before full batch, see README)")
    p.add_argument("--min-pident", type=float, default=85.0,
                    help="min %% identity for the trnH-GUG anchor hit (a tRNA gene -- highly "
                         "conserved across taxa, so a strict default is appropriate)")
    p.add_argument("--min-len", type=int, default=20, help="min alignment length for the trnH-GUG anchor hit")
    p.add_argument("--ndhf-min-pident", type=float, default=70.0,
                    help="min %% identity for the ndhF anchor hit (a protein-coding gene -- "
                         "diverges faster than trnH-GUG across taxa, so this is deliberately "
                         "more permissive; real ndhF orthologs commonly come in around 80-85%% "
                         "nucleotide identity even at moderate taxonomic distance)")
    p.add_argument("--ndhf-min-len", type=int, default=300,
                    help="min alignment length for the ndhF anchor hit (kept well above the "
                         "trnH minimum, since a lower identity threshold needs a longer match "
                         "to stay a confident call rather than noise)")
    p.add_argument("--min-ir-len", type=int, default=1000, help="min alignment length to call a self-BLAST hit an IR copy")
    p.add_argument("--ir-min-pident", type=float, default=95.0, help="min %% identity to call a self-BLAST hit an IR copy")
    p.add_argument("--keep-tmp", action="store_true", help="keep per-sample BLAST scratch files (for debugging)")
    args = p.parse_args()

    if shutil.which("blastn") is None:
        sys.exit("ERROR: blastn not found on PATH. Load/activate a conda env with NCBI BLAST+ first.")

    os.makedirs(args.outdir, exist_ok=True)
    inputs = collect_inputs(args)
    if not inputs:
        sys.exit("No input files found.")

    tmp_dir = tempfile.mkdtemp(prefix="rotate_normalize_")
    log_path = os.path.join(args.outdir, "rotate_normalize_log.tsv")
    fieldnames = [
        "sample", "status", "length", "trnH_anchor_pos_before_rotation", "trnH_orig_strand",
        "ir_len", "ir_pident", "ir_diag_best_len", "ir_diag_best_pident", "ssc_start",
        "ssc_end", "ssc_len", "ndhF_strand_observed", "ssc_flipped", "ndhF_diag_best_pident",
        "ndhF_diag_best_pos", "ndhF_diag_best_strand", "ndhF_diag_in_ssc", "warnings", "output",
    ]

    n_ok, n_fail, n_warn = 0, 0, 0
    with open(log_path, "w", newline="") as log_fh:
        writer = csv.DictWriter(log_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for path in inputs:
            sample = os.path.splitext(os.path.basename(path))[0]
            try:
                row = process_one(path, args, tmp_dir)
                if row["warnings"]:
                    n_warn += 1
                    print(f"[WARN] {sample}: {row['warnings']}", file=sys.stderr)
                else:
                    print(f"[OK]   {sample}: rotated + SSC-normalized ({row['length']} bp)")
                n_ok += 1
            except SampleError as e:
                row = {fn: "" for fn in fieldnames}
                row.update({"sample": sample, "status": f"FAILED: {e}"})
                n_fail += 1
                print(f"[FAIL] {sample}: {e}", file=sys.stderr)
            writer.writerow(row)

    if not args.keep_tmp:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"\nDone: {n_ok} processed ({n_warn} with warnings), {n_fail} failed. Log: {log_path}")


if __name__ == "__main__":
    main()
