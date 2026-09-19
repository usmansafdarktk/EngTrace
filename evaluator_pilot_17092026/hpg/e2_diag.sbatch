#!/bin/bash
# E2 diagnostic: why does Qwen2.5-Math-PRM-72B miss its card's example rewards?
#
#   sbatch -J e2_diag --gpus=2 e2_diag.sbatch
#
# Smoke runs 3 and 4 gave BIT-IDENTICAL self-check rewards under transformers 4.57.3 and
# 4.46.3 ([0.993, 0.0037, 0.159, 0.561] vs the card's [0.992, 0.0048, 0.324, 0.820]),
# with the rotary buffer fp32 in the second - so the library version is not the cause.
# This job scores ONLY the card example (no traces), under each condition that could:
#   qwen7 (default attention)  its card prints rewards too; if it matches, the prompt,
#                              step token and extraction are right and the 72B gap is
#                              numerical, not a pipeline error
#   qwen72 sdpa / eager / flash_attention_2 (if installed)
# Every check also records the rewards with the softmax taken in bf16, as the card does.
# Same scratch rules as e2.sbatch: weights on node-local /tmp, trap-deleted at exit.
#SBATCH --account=fire-finai
#SBATCH --qos=fire-finai
#SBATCH --partition=hpg-rtx6000
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --time=01:00:00
#SBATCH --output=/blue/fire-finai/zhuohan.xie/engtrace-e2/logs/%x_%j.out
set -uo pipefail
BASE=/blue/fire-finai/zhuohan.xie/engtrace-e2
SIF=/blue/fire-finai/zhuohan.xie/finpersona/containers/vllm-openai.sif
SCR=${TMPDIR:-/tmp}/engtrace-e2-$SLURM_JOB_ID
mkdir -p $SCR/hf $SCR/pylib
trap 'rm -rf "$SCR"; echo "cleaned $SCR"' EXIT
echo "job $SLURM_JOB_ID on $(hostname), start $(date -u +%FT%TZ)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
apptainer exec $SIF python3 -m pip install -q --no-deps --target $SCR/pylib \
    transformers==4.46.3 tokenizers==0.20.3 huggingface_hub==0.36.0 2>&1 | tail -3
RUN="apptainer exec --nv --bind $BASE:$BASE --env PYTHONPATH=$SCR/pylib --env HF_HOME=$SCR/hf \
     --env HF_HUB_DISABLE_TELEMETRY=1 --env PYTHONUNBUFFERED=1 $SIF"
echo "flash_attn installed?"; $RUN python3 -c "import flash_attn; print(flash_attn.__version__)" 2>&1 | tail -1

diag() {  # name prm [attn]
    echo "=== $1"
    $RUN python3 $BASE/e2_score.py --prm $2 --inputs $BASE/e2_inputs.jsonl \
        --out $BASE/out/diag_$1.jsonl --selfcheck-only ${3:+--attn $3} 2>&1 | grep -E "SELFCHECK|Error|error:"
}
diag qwen7 qwen7
rm -rf $SCR/hf/hub/models--Qwen--Qwen2.5-Math-PRM-7B
diag qwen72_sdpa qwen72 sdpa
diag qwen72_eager qwen72 eager
diag qwen72_fa2 qwen72 flash_attention_2
echo "end $(date -u +%FT%TZ)"
