#!/bin/bash
# E2 on HiPerGator: one PRM over the E2 inputs, weights on node-local scratch.
#
#   sbatch -J e2_1_smoke --gpus=1 e2.sbatch versa  5    # smoke: self-check + 5 traces
#   sbatch -J e2_2_full  --gpus=2 e2.sbatch qwen72 all  # full run
#
# Job names carry no model names (shared account): 1 = versa, 2 = qwen72, 3 = qwen7.
#
# /blue/fire-finai was 99% full when this was written (163 GB free of 10 TB), so
# NOTHING large goes there: model weights and the few pinned packages live in the
# job's node-local $TMPDIR, which Slurm deletes when the job ends. Only the inputs,
# this script, the scorer and the result files (a few MB) sit in $BASE.
#
# No API key, token or credential is used or placed anywhere: every model is public
# and ungated, and inference is local.
#SBATCH --account=fire-finai
#SBATCH --qos=fire-finai
#SBATCH --partition=hpg-rtx6000
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=03:00:00
#SBATCH --output=/blue/fire-finai/zhuohan.xie/engtrace-e2/logs/%x_%j.out
set -uo pipefail
BASE=/blue/fire-finai/zhuohan.xie/engtrace-e2
SIF=/blue/fire-finai/zhuohan.xie/finpersona/containers/vllm-openai.sif
PRM=$1
LIMIT=${2:-5}

# $TMPDIR on these nodes is the node's SHARED /tmp (3.5 TB local NVMe), not a per-job
# directory Slurm cleans up. So the job makes its own subdirectory and a trap deletes
# it however the job ends - otherwise 146 GB of 72B weights would outlive the job on
# a node other people use.
SCR=${TMPDIR:-/tmp}/engtrace-e2-$SLURM_JOB_ID
mkdir -p $SCR/hf $SCR/pylib
trap 'rm -rf "$SCR"; echo "cleaned $SCR"' EXIT
echo "job $SLURM_JOB_ID on $(hostname), prm=$PRM limit=$LIMIT, start $(date -u +%FT%TZ)"
echo "node-local scratch $SCR:"; df -h $SCR | tail -1
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# Pinned packages, per PRM, installed --no-deps so pip cannot pull in a second torch
# that would shadow the container's Blackwell-capable one.
#  qwen*  transformers 4.46.3 - from the months the PRM was released (its config says
#         4.37, its card >=4.40). Its own modeling_qwen2_rm.py builds rotary positions
#         as arange().type_as(inv_freq); under 4.57.3 that buffer came out bf16, and the
#         self-check failed with an error growing along the sequence (0.12 -> 1.28 in
#         logit over the card's 4 steps). 4.46.3 also still has the cache API it calls.
#  versa  transformers 4.57.3 (E0's pin; VersaPRM runs on transformers' own Llama code,
#         not on custom code) plus peft for its adapter. The container's tokenizers
#         0.23.2 is outside 4.57.3's accepted range, so it is pinned too.
case $PRM in
  qwen*) PKGS="transformers==4.46.3 tokenizers==0.20.3 huggingface_hub==0.36.0" ;;
  *)     PKGS="transformers==4.57.3 tokenizers==0.22.2 huggingface_hub==0.36.0 peft==0.17.1" ;;
esac
echo "packages: $PKGS"
apptainer exec $SIF python3 -m pip install -q --no-deps --target $SCR/pylib $PKGS 2>&1 | tail -3

# Qwen runs with eager attention: of the kernels available here it lands closest to the
# 72B card's example rewards (max gap 0.21 vs 0.26 under sdpa; hpg/README.md).
case $PRM in qwen*) ATTN="--attn eager" ;; *) ATTN="" ;; esac
if [ "$LIMIT" = "all" ]; then LIM=""; TAG=full; else LIM="--limit $LIMIT"; TAG=smoke$LIMIT; fi
OUT=$BASE/out/${PRM}_${TAG}.jsonl
# No explicit bind for $SCR: Apptainer mounts the host's /tmp by default (the pip
# step above writes there with no bind at all). Binding a SUBDIRECTORY of /tmp fails
# on this image - it ships its own /tmp/ep-kernels-wheels.sha256, and the underlay
# mount of the new path collides with it (FATAL: mount hook function failure).
apptainer exec --nv --bind $BASE:$BASE \
    --env PYTHONPATH=$SCR/pylib --env HF_HOME=$SCR/hf --env HF_HUB_DISABLE_TELEMETRY=1 \
    --env PYTHONUNBUFFERED=1 $SIF \
    python3 $BASE/e2_score.py --prm $PRM --inputs $BASE/e2_inputs.jsonl --out $OUT $LIM $ATTN
RC=$?
echo "rc=$RC end=$(date -u +%FT%TZ)"
echo "scratch used:"; du -sh $SCR/hf 2>/dev/null
exit $RC
