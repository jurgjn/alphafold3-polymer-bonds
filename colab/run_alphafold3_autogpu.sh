#!/usr/bin/env bash
# https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md

cd /app/alphafold
git checkout src/alphafold3/model/model_config.py

ARGS_XTRA=""

echo Using compilation time workaround with XLA flags
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"

GPU_CAPABILITY=`nvidia-smi --query-gpu=compute_cap --format=csv,noheader | cut -d '.' -f 1`
if [[ "$GPU_CAPABILITY" == 7 ]] # e.g. T4 GPU is capability 7
then
    echo GPU capability 7, adjusting XLA_FLAGS and --flash-attention_implementation
    export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
    ARGS_XTRA="${ARGS_XTRA} --flash_attention_implementation=xla"
fi

GPU_MEMAVAIL=`nvidia-smi --query-gpu=memory.total --format=csv,noheader | cut -d ' ' -f 1`
if [ "$GPU_MEMAVAIL" -gt 80000 ]
then
    echo "Using default GPU memory setup (80 GB)"
    export XLA_PYTHON_CLIENT_PREALLOCATE=true
    export XLA_CLIENT_MEM_FRACTION=0.95
else
    echo "Using low GPU memory setup (40 GB)"
    echo Adjusting pair_transition_shard_spec in model_config.py
    git apply <<EOF
diff --git a/src/alphafold3/model/model_config.py b/src/alphafold3/model/model_config.py
index 2040d8f..54d13fc 100644
--- a/src/alphafold3/model/model_config.py
+++ b/src/alphafold3/model/model_config.py
@@ -26,7 +26,8 @@ class GlobalConfig(base_config.BaseConfig):
   pair_attention_chunk_size: Sequence[_Shape2DType] = ((1536, 128), (None, 32))
   pair_transition_shard_spec: Sequence[_Shape2DType] = (
       (2048, None),
-      (None, 1024),
+      (3072, 1024),
+      (None, 512),
   )
   # Note: flash_attention_implementation = 'xla' means no flash attention.
   flash_attention_implementation: attention.Implementation = 'triton'
EOF

    echo Enabling unified memory
    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export TF_FORCE_UNIFIED_MEMORY=true
    export XLA_CLIENT_MEM_FRACTION=3.2
fi

echo Executing run_alphafold.py
python3 run_alphafold.py $@ $ARGS_XTRA
