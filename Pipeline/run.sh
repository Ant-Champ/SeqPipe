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
# Without SLURM (a workstation): SEQPIPE_NO_SLURM=1 ./run.sh CONFIG  (uses every CPU, 90 % of free memory)
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${SEQPIPE_ENV:-seqpipe}"                                  # driver env (envs/snakemake.yaml)
CONDA_ENVS_DIR="${SEQPIPE_CONDA_PREFIX:-$PIPELINE_DIR/.conda_envs}" # per-rule envs, shared by all runs

usage() { sed -n '2,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1; }
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
    # setup: the depth check (a checkpoint) is switched off so that snakemake sees
    # every step of the test config -- otherwise only the steps before it are
    # known and their conda environments alone would be built
    setup) config="$PIPELINE_DIR/test/config.yaml"; extra=(--conda-create-envs-only --config min_theoretical_depth=0) ;;
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
# On a workstation (no SLURM), set SEQPIPE_NO_SLURM=1 to use this machine's CPUs and memory.
if [ "$mode" != setup ] && [ "$dry" = no ] && [ -z "${SLURM_JOB_ID:-}" ] && [ "${SEQPIPE_NO_SLURM:-0}" != 1 ]; then
    cat >&2 <<'EOF'
ERROR: not inside a SLURM allocation (never run analyses on a login node).
Start one first, from inside tmux so it survives a dropped connection, e.g.:
    tmux new -s pipeline
    salloc --account=<your-account> --time=4:00:00 --cpus-per-task=16 --mem=64G
then re-run this command. (Add -n for a dry run, which is fine on a login node.)
Not on a cluster? Run it as:  SEQPIPE_NO_SLURM=1 ./run.sh CONFIG
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

# --- every conda environment this run needs must already be built ------------
# Compute nodes have no internet: an environment missing there makes conda wait
# forever (the run "stalls" when it reaches that step). Check before starting.
if [ "$mode" != setup ] && [ "$dry" = no ]; then
    # (a failure here must be shown, not end the script silently: set -e / pipefail)
    if ! envlist=$(snakemake --snakefile "$PIPELINE_DIR/Snakefile" --configfile "$config" --directory "$run_dir" \
                       --use-conda --conda-prefix "$CONDA_ENVS_DIR" --cores 1 --list-conda-envs \
                       --config min_theoretical_depth=0 2>&1); then
        echo "ERROR: snakemake could not read the workflow (config, sample sheet, references...):" >&2
        printf '%s\n' "$envlist" | grep -v "^Building DAG\|^Assuming unrestricted\|^host:" | tail -n 40 >&2
        exit 1
    fi
    missing=$(printf '%s\n' "$envlist" \
              | awk -F'\t' 'NF >= 3 && $1 != "environment" && $3 ~ /\// { print $1 "\t" $3 }' \
              | while IFS=$'\t' read -r env loc; do [ -d "$loc" ] || echo "  $(basename "$env")"; done || true)
    if [ -n "$missing" ]; then
        cat >&2 <<EOF
ERROR: these conda environments are not built yet:
$missing
Compute nodes have no internet, so the run would hang when it reaches them.
Build them once on a login node (after every 'git pull' that changes envs/):
    $PIPELINE_DIR/run.sh setup
EOF
        exit 1
    fi
fi

echo "+ ${cmd[*]}"
"${cmd[@]}"

if [ "$mode" = setup ]; then
    echo "All conda environments are built in $CONDA_ENVS_DIR"
elif [ "$dry" = no ]; then
    out=$(grep -E '^out_dir' "$config" | sed -E 's/^out_dir *: *"?([^"]*)"?.*/\1/')
    echo "Done. Start with: ${out}report.html (copy it to your PC and open it in a browser) or ${out}summary.tsv"
fi
