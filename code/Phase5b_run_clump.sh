#!/bin/bash
# 06b_run_clump.sh — 对每个基因的 SNP 列表运行 PLINK LD clump
# 用法: bash scripts/06b_run_clump.sh

set -e

PLINK="<path_to_1000G_reference>/tools/plink.exe"
BFILE="<path_to_1000G_reference>/data/ref/1000G_EUR/EUR"
INPUT_DIR="data/clump_inputs"
OUTPUT_DIR="data/clump_outputs"

mkdir -p "$OUTPUT_DIR"

# 清理旧输出
rm -f "$OUTPUT_DIR"/*.clumped "$OUTPUT_DIR"/*.log "$OUTPUT_DIR"/*.nosex 2>/dev/null || true

count=0
total=$(ls "$INPUT_DIR"/*.txt 2>/dev/null | wc -l)
echo "=== PLINK LD Clump: $total 个基因 ==="

for input_file in "$INPUT_DIR"/*.txt; do
  gene=$(basename "$input_file" .txt)
  output_prefix="$OUTPUT_DIR/$gene"

  "$PLINK" \
    --bfile "$BFILE" \
    --clump "$input_file" \
    --clump-kb 10000 \
    --clump-r2 0.001 \
    --clump-p1 1 \
    --clump-p2 1 \
    --out "$output_prefix" \
    > /dev/null 2>&1 || true

  count=$((count + 1))
  if [ $((count % 20)) -eq 0 ]; then
    echo "  进度: $count / $total"
  fi
done

echo "=== Clump 完成: $count 个基因已处理 ==="
echo "成功生成 clumped 文件数: $(ls "$OUTPUT_DIR"/*.clumped 2>/dev/null | wc -l)"
