#!/bin/bash
#
# mitobim_rescue.sh
#
# MITObim (quick mode) rescue assembly for a sample whose GetOrganelle (or
# other) assembly came out fragmented rather than a complete, circular
# plastome. Seeds MITObim with a reference plastome from a related taxon
# and extends using the sample's own trimmed reads plus its existing
# fragmented scaffold(s), for up to 10 rounds.
#
# Requires: MIRA (v4.0.2 used in this project; HPC module or on PATH) and
# MITObim.pl (v1.9.1; typically installed in its own conda env, since it
# depends on an older Perl/MIRA toolchain).
#
# ------------------------------------------------------------------------
# USAGE EXAMPLE (single sample)
# ------------------------------------------------------------------------
#   ./mitobim_rescue.sh \
#       --sample SRR12345678 \
#       --seed /path/to/ref_cps/Calendula_arvensis_cp.fasta \
#       --r1 SRR12345678_1.trimmed.fastq \
#       --r2 SRR12345678_2.trimmed.fastq \
#       --scaffold SRR12345678.path_sequence.fasta \
#       --outdir mitobim_out/SRR12345678 \
#       --mira-module mira/4.0.2 \
#       --conda-env mitobim
#
# --scaffold can be repeated for samples with more than one GetOrganelle
# scaffold. To process a LIST of samples, loop over this script, e.g.:
#
#   while read -r sample ref r1 r2 scaffold; do
#       ./mitobim_rescue.sh --sample "$sample" --seed "$ref" \
#           --r1 "$r1" --r2 "$r2" --scaffold "$scaffold" \
#           --outdir "mitobim_out/${sample}"
#   done < my_sample_list.tsv
#
# ------------------------------------------------------------------------

set -euo pipefail

SAMPLE=""
SEED_FASTA=""
R1=""
R2=""
SCAFFOLDS=()
OUTDIR=""
MIRA_MODULE="mira/4.0.2"
CONDA_ENV_NAME="mitobim"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift 2 ;;
        --seed) SEED_FASTA="$2"; shift 2 ;;
        --r1) R1="$2"; shift 2 ;;
        --r2) R2="$2"; shift 2 ;;
        --scaffold) SCAFFOLDS+=("$2"); shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --mira-module) MIRA_MODULE="$2"; shift 2 ;;
        --conda-env) CONDA_ENV_NAME="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

for req_name in SAMPLE SEED_FASTA R1 R2 OUTDIR; do
    if [[ -z "${!req_name}" ]]; then
        echo "ERROR: --${req_name,,} is required" >&2
        exit 1
    fi
done
if [[ ${#SCAFFOLDS[@]} -eq 0 ]]; then
    echo "ERROR: at least one --scaffold is required" >&2
    exit 1
fi

# mirabait (part of MIRA) hits a glibc locale assertion and crashes
# silently on some HPC nodes ("loadlocale.c: _nl_intern_locale_data:
# Assertion ... failed"), which MITObim then misreports as "no reads
# match your reference". Forcing a plain, always-available locale avoids
# it (see https://github.com/chrishah/MITObim/issues/49).
export LC_ALL=C
export LANG=C

mkdir -p "${OUTDIR}"

ln -sf "$(readlink -f "${R1}")" "${OUTDIR}/${SAMPLE}_R1_001_P.fastq"
ln -sf "$(readlink -f "${R2}")" "${OUTDIR}/${SAMPLE}_R2_001_P.fastq"

scaffold_names=()
for f in "${SCAFFOLDS[@]}"; do
    bn=$(basename "${f}")
    ln -sf "$(readlink -f "${f}")" "${OUTDIR}/${bn}"
    scaffold_names+=("${bn}")
done

cp -f "${SEED_FASTA}" "${OUTDIR}/seedOrganelle.fasta"

# ---- MIRA (HPC module) + MITObim (conda env) -------------------------------
# Lmod's own init script (and some conda init scripts) reference variables
# like $LD_LIBRARY_PATH without a default and aren't `set -u` safe, so relax
# nounset just for this block and restore it right after.
set +u
if ! command -v module >/dev/null 2>&1; then
    for f in /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
        [[ -f "${f}" ]] && source "${f}" && break
    done
fi
if command -v module >/dev/null 2>&1; then
    module load "${MIRA_MODULE}"
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
set -u

mira_bin=$(command -v mira || true)
if [[ -z "${mira_bin}" ]]; then
    echo "ERROR: 'mira' not found on PATH after 'module load ${MIRA_MODULE}'" >&2
    exit 1
fi
mira_dir=$(dirname "${mira_bin}")
echo "mirapath: ${mira_dir}"

mitobim_bin=$(command -v MITObim.pl || true)
if [[ -z "${mitobim_bin}" ]]; then
    echo "ERROR: 'MITObim.pl' not found on PATH after activating conda env '${CONDA_ENV_NAME}'" >&2
    exit 1
fi

# ---- run MITObim, quick mode, 10 rounds max --------------------------------
cd "${OUTDIR}"

"${mitobim_bin}" \
    -start 1 -end 10 \
    -sample "${SAMPLE}" \
    -ref "${SAMPLE}" \
    --quick seedOrganelle.fasta \
    -readpool "${SAMPLE}_R1_001_P.fastq" "${SAMPLE}_R2_001_P.fastq" "${scaffold_names[@]}" \
    -mirapath "${mira_dir}" \
    --clean --verbose \
    > "${OUTDIR}/${SAMPLE}_mitobim.log" 2>&1

status=$?
echo "MITObim.pl exited with status ${status} for ${SAMPLE} (see ${OUTDIR}/${SAMPLE}_mitobim.log)"
exit ${status}


