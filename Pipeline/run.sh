#!/usr/bin/env bash
# Launcher for the pipeline. Usage:
#
#   ./run.sh setup                          build every conda environment (login node, once)
#   ./run.sh test   [snakemake options]     run the bundled test data (inside salloc)
#   ./run.sh CONFIG [snakemake options]     run an analysis (inside salloc)
#
# Useful snakemake options: -n (dry run: list jobs, run nothing; OK on a login node),
#   --forcerun lowcov (redo that step), --until align_sr (stop after it), --unlock.
# Uses every CPU and all memory of the current salloc allocation.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${SEQPIPE_ENV:-seqpipe}"                                  # driver env (envs/snakemake.yaml)
CONDA_ENVS_DIR="${SEQPIPE_CONDA_PREFIX:-$PIPELINE_DIR/.conda_envs}" # per-rule envs, shared by all runs

usage() { sed -n '2,10p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1; }
[ $# -ge 1 ] || usage

# --- driver environment (snakemake) ----------------------------------------
if [ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]; then
    if ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: conda not found. Load/initialise conda first (see README, Installation)." >&2
        exit 1
    fi
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME" || {
        echo "ERROR: conda environment '$ENV_NAME' not found. Create it once with:" >&2
        echo "  conda env create -f $PIPELINE_DIR/envs/snakemake.yaml" >&2
        exit 1
    }
fi

# --- what to run -------------------------------------------------------------
mode="$1"; shift
extra=()
case "$mode" in
    -h|--help) usage ;;
    setup) config="$PIPELINE_DIR/test/config.yaml"; extra=(--conda-create-envs-only) ;;
    test)  config="$PIPELINE_DIR/test/config.yaml" ;;
    *)     config="$(readlink -f "$mode")" ;;
esac
[ -f "$config" ] || { echo "ERROR: config file not found: $mode" >&2; exit 1; }
run_dir="$(dirname "$config")"     # snakemake bookkeeping (.snakemake/) and relative paths live here

dry=no
for a in "$@"; do
    case "$a" in -n|--dry-run|--dryrun) dry=yes ;; esac
done

# --- resources: the salloc allocation ------------------------------------------
if [ "$mode" != setup ] && [ "$dry" = no ] && [ -z "${SLURM_JOB_ID:-}" ]; then
    cat >&2 <<'EOF'
ERROR: not inside a SLURM allocation (never run analyses on a login node).
Start one first, from inside tmux so it survives a dropped connection, e.g.:
    tmux new -s pipeline
    salloc --account=def-rodrigu1 --time=4:00:00 --cpus-per-task=16 --mem=64G
then re-run this command. (Add -n for a dry run, which is fine on a login node.)
EOF
    exit 1
fi

if [ "$mode" = setup ]; then
    cores=1
else
    cores="${SLURM_CPUS_PER_TASK:-$(nproc)}"
fi
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    mem_mb="$SLURM_MEM_PER_NODE"
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    mem_mb=$(( SLURM_MEM_PER_CPU * cores ))
else
    mem_mb=$(awk '/MemAvailable/ {print int($2 / 1024)}' /proc/meminfo)
fi
mem_mb=$(( mem_mb * 9 / 10 ))      # headroom for snakemake itself
export TMPDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"

cmd=(snakemake
     --snakefile "$PIPELINE_DIR/Snakefile"
     --configfile "$config"
     --directory "$run_dir"
     --use-conda --conda-prefix "$CONDA_ENVS_DIR"
     --cores "$cores" --resources mem_mb="$mem_mb"
     --rerun-incomplete --keep-going --printshellcmds
     ${extra[@]+"${extra[@]}"} "$@")

echo "+ ${cmd[*]}"
"${cmd[@]}"

if [ "$mode" = setup ]; then
    echo "All conda environments are built in $CONDA_ENVS_DIR"
elif [ "$dry" = no ]; then
    echo "Done. Start with: $(grep -E '^out_dir' "$config" | sed -E 's/^out_dir *: *"?([^"]*)"?.*/\1/')summary.tsv"
fi
