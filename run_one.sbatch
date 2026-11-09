#!/usr/bin/env bash
#SBATCH --job-name=moaz04S1_ica
#SBATCH --account=def-jmeltzer
#SBATCH --output=/home/vant7ee/scratch/FacerecognitionData/logs/moaz04S1_%j.out
#SBATCH --error=/home/vant7ee/scratch/FacerecognitionData/logs/moaz04S1_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

set -euo pipefail
module load "${MATLAB_MODULE:-matlab}"
BASE="/home/vant7ee/scratch/FacerecognitionData"
# Slurm runs a spool copy: do not use dirname(BASH_SOURCE) here.
export MOAZ_SCRIPT_DIR="${MOAZ_SCRIPT_DIR:-$BASE/preprocessing_moaz_native}"
export MOAZ_LW_ROOT="${MOAZ_LW_ROOT:-$BASE/letswave6_moaz_pin}"
export MOAZ_BDF="${MOAZ_BDF:-$BASE/bdf/subj04_S1.bdf}"
MOAZ_RUN_ID="${SLURM_JOB_ID:-interactive_$(date +%Y%m%d_%H%M%S)_$$}"
export MOAZ_OUTPUT="${MOAZ_OUTPUT:-$BASE/preprocessed_moaz_native_${MOAZ_RUN_ID}/subj04/S1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
MOAZ_TEMP_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
MOAZ_CONFIG_DIR="$(mktemp -d "$MOAZ_TEMP_BASE/moaz_config.XXXXXX")"
export MOAZ_CONFIG="$MOAZ_CONFIG_DIR/moaz_config.json"
trap 'rm -f -- "$MOAZ_CONFIG"; rmdir -- "$MOAZ_CONFIG_DIR"' EXIT
test -f "$MOAZ_SCRIPT_DIR/moaz_diagnose.m"
test -d "$MOAZ_LW_ROOT"
test -f "$MOAZ_BDF"
# Environment variables are read inside Matlab, never interpolated as Matlab code.
matlab -batch "addpath(getenv('MOAZ_SCRIPT_DIR')); moaz_entrypoint;"
