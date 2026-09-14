#!/usr/bin/env python3
"""
novoplasty_to_gfa.py

Builds Bandage-readable GFA (v1) files from NOVOPlasty output.

For each sample it uses:
  - Contigs_1_<sample>.fasta            -> GFA segments (S lines / nodes)
  - Merged_contigs_<sample>.txt         -> parsed for the "LINKS BETWEEN CONTIGS" section
  - Circularized_assembly_1_<sample>.fasta -> if present, signals that a single-contig assembly was confirmed circular; used to add a self-loop link so Bandage draws the contig as a circle

Usage:
------------------------------------------------------------------------
Single sample:
    python novoplasty_to_gfa.py Contigs_1_ERR10116773.fasta Merged_contigs_ERR10116773.txt out.gfa
or
    python novoplasty_to_gfa.py Circularized_assembly_1_ERR13996271.fasta out.gfa

Batch, e.g. over a directory of per-sample NOVOPlasty output folders:
    python novoplasty_to_gfa.py --batch /path/to/novoplasty_results

    Recursively finds every Contigs_1_*.fasta under the root, pairs it with whatever Merged_contigs_*.txt / Circularized_assembly_1_*.fasta exist in the same folder (both optional), and writes <sample>.gfa next to the inputs. Also separately finds every Circularized_assembly_1_*.fasta that has no matching Contigs_1 file (first-pass circularizations) and converts those too, using the circularized fasta itself as the contig source. A summary TSV is written at the end (default: novoplasty_to_gfa_summary.tsv in the current directory).

    Options:
      --out-dir DIR   write all .gfa files here instead of next to inputs
      --summary FILE  path for the summary TSV (default ./novoplasty_to_gfa_summary.tsv)
------------------------------------------------------------------------
"""

import argparse
import re
import sys
from pathlib import Path


def parse_fasta(path):
    """Return dict: contig_number(str, normalized no leading zeros) -> (header, sequence)."""
    contigs = {}
    header = None
    seq_chunks = []

    def flush():
        if header is not None:
            seq = "".join(seq_chunks)
            num = extract_contig_number(header)
            key = num if num is not None else header
            contigs[key] = (header, seq)

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        flush()

    return contigs


def read_first_fasta_seq(path):
    """Return just the sequence of the first record in a fasta file (used for Circularized_assembly files)."""
    seq_chunks = []
    started = False
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            seq_chunks.append(line.strip())
    return "".join(seq_chunks)


def extract_contig_number(text):
    """Pull the numeric contig id out of a header or log line, stripping leading zeros."""
    m = re.search(r"[Cc]ontig[_\s]?0*(\d+)", text)
    if m:
        return str(int(m.group(1)))
    return None


def find_circularized_fasta(folder, sample):
    candidates = sorted(Path(folder).glob(f"Circularized_assembly*{sample}*.fasta"))
    return candidates[0] if candidates else None


def _norm_node(token):
    token = token.strip()
    if token.upper() in ("START", "END"):
        return token.upper()
    m = re.match(r"0*(\d+)$", token)
    return m.group(1) if m else token


def parse_merged_contigs(path):
    """
    Parse the real "LINKS BETWEEN CONTIGS" + "OPTION N" format.

    Returns:
      options: list of dicts {option, arrangement (list of contig ids), length}
      links:   list of dicts {a, b, overlap_bp (int or None), evidence (str)}
               -- includes both direct links from the LINKS block and
               inferred circular-closure links (END contigs back to
               START contigs), the latter always with overlap_bp=None
               so they get flagged NEEDS REVIEW downstream.
    """
    text = Path(path).read_text(errors="replace")
    # NOTE: don't truncate at the first fasta record -- OPTION 2, 3, ... blocks
    # each come after their own full-genome sequence, so later options would
    # be silently skipped if we cut the text short here.

    # (\S+?) is lazy so it stops BEFORE the dashes instead of swallowing them --
    # with a greedy \S+ here, "01----> 02" mis-captures the source as "01---"
    # (extra dashes attached) because \S+ backtracks the minimum amount needed
    # to let -+> match, not the amount that gives the "right" token boundary.
    edge_re = re.compile(r"^(\S+?)-+>\s*(.+?)\s*$", re.MULTILINE)
    raw_edges = []
    for m in edge_re.finditer(text):
        src = _norm_node(m.group(1))
        dsts = [_norm_node(d) for d in re.split(r"\s+OR\s+", m.group(2))]
        raw_edges.append((src, dsts))

    start_targets = set()
    end_sources = set()
    links = []
    seen = set()

    for src, dsts in raw_edges:
        for dst in dsts:
            if src == "START":
                start_targets.add(dst)
                continue
            if dst == "END":
                end_sources.add(src)
                continue
            key = tuple(sorted((src, dst)))
            if key in seen:
                continue
            seen.add(key)
            links.append({
                "a": src, "b": dst, "overlap_bp": None,
                "evidence": f"LINKS BETWEEN CONTIGS: {src}----> {dst}",
            })

    for s in end_sources:
        for t in start_targets:
            closure_key = (s, t, "closure")
            if closure_key in seen:
                continue
            seen.add(closure_key)
            links.append({
                "a": s, "b": t, "overlap_bp": None,
                "evidence": f"inferred circular closure ({s}----> END, START---->{t}); NOVOPlasty did not state this explicitly",
            })

    options = []
    for om in re.finditer(
        r"OPTION\s+(\d+).*?Contig Arrangement\s*:\s*(.+?)\s*\n\s*Assembly length\s*:\s*(\d+)\s*bp",
        text, re.DOTALL,
    ):
        opt_num, arrangement, length = om.groups()
        contigs_in_opt = [_norm_node(c) for c in re.split(r"\+", arrangement.strip())]
        options.append({"option": opt_num, "arrangement": contigs_in_opt, "length": int(length)})

    return options, links


def attach_overlaps_from_options(contigs, options, links):
    """
    For links between exactly two contigs, if some OPTION's arrangement is
    exactly that same pair, back-calculate overlap_bp = len(A)+len(B) -
    option_length and attach it (upgrading the link from an unknown/0M
    CIGAR to a real one). Sanity-bounded the same way as the circularized
    self-loop calculation.
    """
    pair_lengths = {}
    for opt in options:
        if len(opt["arrangement"]) == 2:
            key = tuple(sorted(opt["arrangement"]))
            pair_lengths[key] = opt["length"]

    for link in links:
        if link["overlap_bp"] is not None:
            continue
        a, b = link["a"], link["b"]
        if a in ("START", "END") or b in ("START", "END"):
            continue
        key = tuple(sorted((a, b)))
        if key not in pair_lengths or a not in contigs or b not in contigs:
            continue
        len_a = len(contigs[a][1])
        len_b = len(contigs[b][1])
        option_len = pair_lengths[key]
        overlap = len_a + len_b - option_len
        sane = 0 < overlap < max(200, min(len_a, len_b) * 0.2)
        if sane:
            link["overlap_bp"] = overlap
            link["evidence"] += f" | overlap inferred from OPTION arrangement {a}+{b} = {option_len} bp"


def add_self_loop_if_circularized(contigs, circ_path, links):
    """If there's exactly one contig and a Circularized_assembly file exists next to it,
    add a self-loop link so Bandage draws the node as a circle. Returns True if a
    self-loop was added."""
    if len(contigs) != 1 or circ_path is None:
        return False

    key = next(iter(contigs))
    orig_seq = contigs[key][1]
    circ_seq = read_first_fasta_seq(circ_path)

    if orig_seq == circ_seq:
        # NOVOPlasty circularized on the first pass and never wrote a separate
        # Contigs_1_<sample>.fasta -- the "contig" we parsed IS the
        # Circularized_assembly file itself, so there's no pre-circularization
        # sequence left to diff against and the overlap size can't be computed.
        # Circularization is still CONFIRMED (NOVOPlasty's own log says so),
        # just flagged for review because the overlap length is unknown.
        links.append({
            "a": key,
            "b": key,
            "overlap_bp": None,
            "evidence": (
                f"self-loop: {circ_path.name} IS the contig source -- NOVOPlasty "
                f"closed the genome circular on the first pass and wrote no "
                f"separate Contigs_1 file, so overlap length can't be back-"
                f"calculated (circularization itself is confirmed, not inferred)"
            ),
        })
        return True

    overlap = len(orig_seq) - len(circ_seq)
    sane = 0 < overlap < max(200, len(orig_seq) * 0.05)

    links.append({
        "a": key,
        "b": key,
        "overlap_bp": overlap if sane else None,
        "evidence": (
            f"self-loop: {circ_path.name} confirms circularization "
            f"(Contigs_1 seq {len(orig_seq)} bp vs circularized seq {len(circ_seq)} bp)"
        ),
    })
    return True


def write_gfa(contigs, links, out_path):
    needs_review = []
    with open(out_path, "w") as out:
        out.write("H\tVN:Z:1.0\n")

        for key, (header, seq) in contigs.items():
            seg_id = f"contig{key}" if key.isdigit() else re.sub(r"\W+", "_", header)
            out.write(f"S\t{seg_id}\t{seq}\tLN:i:{len(seq)}\n")

        for link in links:
            a, b = link["a"], link["b"]
            if a in ("START", "END") or b in ("START", "END"):
                continue
            if a not in contigs or b not in contigs:
                continue
            seg_a = f"contig{a}"
            seg_b = f"contig{b}"
            overlap = link["overlap_bp"]
            cigar = f"{overlap}M" if overlap else "0M"
            out.write(f"L\t{seg_a}\t+\t{seg_b}\t+\t{cigar}\n")
            if overlap is None:
                needs_review.append((seg_a, seg_b, link["evidence"]))

    return needs_review


def convert_one(fasta_path, merged_path, out_path, circ_path=None, verbose=True):
    """
    Run the fasta (+ optional merged log, + optional circularized fasta) -> gfa
    conversion for a single sample. merged_path may be None if that file
    doesn't exist. fasta_path may itself be a Circularized_assembly_1_*.fasta
    (first-pass circularization, no separate Contigs_1 file) -- in that case
    pass circ_path=fasta_path explicitly. Returns a result dict.
    """
    contigs = parse_fasta(fasta_path)

    links = []
    options = []
    if merged_path is not None:
        options, links = parse_merged_contigs(merged_path)
        attach_overlaps_from_options(contigs, options, links)

    if circ_path is None:
        fasta_name = Path(fasta_path).name
        if fasta_name.startswith("Circularized_assembly"):
            # fasta_path IS the circularized fasta (first-pass circularization,
            # no separate Contigs_1 file) -- it's its own circ_path.
            circ_path = Path(fasta_path)
        else:
            circ_path = find_circularized_fasta(Path(fasta_path).parent, Path(fasta_path).stem.replace("Contigs_1_", ""))
    self_loop_added = add_self_loop_if_circularized(contigs, circ_path, links)

    needs_review = write_gfa(contigs, links, out_path)

    if verbose:
        print(f"Parsed {len(contigs)} contigs from {fasta_path}:")
        for key, (header, seq) in contigs.items():
            print(f"  contig{key}\t{len(seq)} bp\t(header: {header})")
        if merged_path is not None:
            print(f"\nOPTION arrangements found in {merged_path}:")
            for opt in options:
                print(f"  Option {opt['option']}: {'+'.join(opt['arrangement'])} = {opt['length']} bp")
            print(f"\nLinks parsed (including inferred circular closures): {len(links)}")
        else:
            print("\nNo Merged_contigs file given/found -- skipping link parsing.")
        if self_loop_added:
            print(f"Circularized_assembly file found ({circ_path.name}) -- added self-loop link.")
        for link in links:
            bp = link["overlap_bp"] if link["overlap_bp"] is not None else "unknown"
            arrow = "(self-loop)" if link["a"] == link["b"] else ""
            print(f"  contig{link['a']} -- contig{link['b']} {arrow}\toverlap={bp}\t| {link['evidence']}")
        print(f"\nWrote {out_path}")
        if needs_review:
            print("\nNEEDS REVIEW -- these links had no explicit/confident overlap length")
            print("(written as 0M; check in Bandage and edit/delete the L line if the connection looks wrong):")
            for seg_a, seg_b, evidence in needs_review:
                print(f"  {seg_a} -- {seg_b}\t| {evidence}")

    return {
        "num_contigs": len(contigs),
        "num_links": len(links),
        "num_needs_review": len(needs_review),
        "self_loop_added": self_loop_added,
    }


def _process_sample(fasta_path, sample, folder, merged_path, circ_path, out_dir, rows):
    dest_dir = Path(out_dir) if out_dir else folder
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_path = dest_dir / f"{sample}.gfa"

    print(f"=== {sample}  ({folder}) ===")

    try:
        result = convert_one(fasta_path, merged_path, out_path, circ_path=circ_path, verbose=False)

        if result["self_loop_added"]:
            status = "ok_circular"
        elif result["num_links"] > 0:
            status = "ok"
        elif result["num_contigs"] == 1:
            status = "ok_single_contig_no_circular_evidence"
        else:
            status = "ok_no_links_found"

        print(
            f"  {result['num_contigs']} contigs, {result['num_links']} links "
            f"({result['num_needs_review']} need review) "
            f"[merged_file={'yes' if merged_path else 'no'}, circularized_file={'yes' if circ_path else 'no'}] "
            f"-> {out_path.name}\n"
        )
        rows.append([
            str(folder), sample,
            result["num_contigs"], result["num_links"], result["num_needs_review"],
            "yes" if merged_path else "no", "yes" if circ_path else "no",
            status, "",
        ])
    except Exception as e:
        print(f"  ERROR: {e}\n")
        rows.append([str(folder), sample, "", "", "", "", "", "error", str(e)])


def run_batch(root, out_dir, summary_path):
    root = Path(root)
    fasta_files = sorted(root.rglob("Contigs_1_*.fasta"))
    contigs_samples = {f.name[len("Contigs_1_"):-len(".fasta")] for f in fasta_files}

    # Samples that circularized on NOVOPlasty's first assembly pass never get
    # a Contigs_1_<sample>.fasta -- only Circularized_assembly_1_<sample>.fasta
    # is written. Find those separately so they aren't silently dropped.
    circ_only = []
    for circ_path in sorted(root.rglob("Circularized_assembly_1_*.fasta")):
        sample = circ_path.name[len("Circularized_assembly_1_"):-len(".fasta")]
        if sample not in contigs_samples:
            circ_only.append((sample, circ_path))

    if not fasta_files and not circ_only:
        print(f"No Contigs_1_*.fasta or Circularized_assembly_1_*.fasta files found under {root}")
        return

    print(
        f"Found {len(fasta_files)} sample(s) with Contigs_1 files and "
        f"{len(circ_only)} sample(s) circularized on the first pass "
        f"(Circularized_assembly only, no Contigs_1) under {root}\n"
    )

    rows = []
    for fasta_path in fasta_files:
        folder = fasta_path.parent
        sample = fasta_path.name[len("Contigs_1_"):-len(".fasta")]
        merged_path = folder / f"Merged_contigs_{sample}.txt"
        merged_path = merged_path if merged_path.exists() else None
        circ_path = find_circularized_fasta(folder, sample)
        _process_sample(fasta_path, sample, folder, merged_path, circ_path, out_dir, rows)

    for sample, circ_path in circ_only:
        folder = circ_path.parent
        merged_path = folder / f"Merged_contigs_{sample}.txt"
        merged_path = merged_path if merged_path.exists() else None
        # circ_path IS the fasta source here -- pass it as both so
        # add_self_loop_if_circularized recognizes orig_seq == circ_seq.
        _process_sample(circ_path, sample, folder, merged_path, circ_path, out_dir, rows)

    summary_path = Path(summary_path)
    with open(summary_path, "w") as f:
        f.write(
            "folder\tsample\tnum_contigs\tnum_links\tnum_needs_review\t"
            "merged_file_found\tcircularized_file_found\tstatus\tmessage\n"
        )
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    def count(status_prefix):
        return sum(1 for r in rows if str(r[7]).startswith(status_prefix))

    print("=== Summary ===")
    print(f"  {count('ok_circular')} confirmed circular (single contig, self-loop added)")
    print(f"  {count('ok') - count('ok_circular') - count('ok_no_links_found') - count('ok_single_contig')} multi-contig with merge links")
    print(f"  {count('ok_no_links_found')} multi-contig, no links detected")
    print(f"  {count('ok_single_contig_no_circular_evidence')} single contig, circularity unconfirmed")
    print(f"  {count('error')} errored")
    print(f"  Details written to {summary_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fasta", nargs="?", help="Contigs_1_*.fasta, or Circularized_assembly_1_*.fasta if no Contigs_1 file exists (single-sample mode)")
    parser.add_argument("merged", nargs="?", help="Merged_contigs_*.txt, or '-' if it doesn't exist (single-sample mode)")
    parser.add_argument("out", nargs="?", help="output .gfa path (single-sample mode)")
    parser.add_argument("--batch", metavar="ROOT_DIR", help="recursively convert every sample under this directory")
    parser.add_argument("--out-dir", metavar="DIR", help="[batch mode] write all .gfa files here instead of next to inputs")
    parser.add_argument("--summary", metavar="FILE", default="novoplasty_to_gfa_summary.tsv",
                         help="[batch mode] path for the summary TSV (default: %(default)s)")
    args = parser.parse_args()

    if args.batch:
        run_batch(args.batch, args.out_dir, args.summary)
    elif args.fasta and args.out:
        merged_path = None if (args.merged in (None, "-")) else args.merged
        convert_one(args.fasta, merged_path, args.out, verbose=True)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
