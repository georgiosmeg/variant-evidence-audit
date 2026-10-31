#!/usr/bin/env bash
set -euo pipefail

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release

V421=HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.clean.bed
V50Q=HG002_GRCh38_v5.0q_smvar.benchmark.clean.bed

for f in "$V421" "$V50Q" clinvar_plp.bed lost_in_v50q.bed gained_in_v50q.bed; do
  [ -s "$f" ] || { echo "missing $f, run 01 first" >&2; exit 1; }
done

echo "=== A. per-gene proportions ==="
printf "gene\ttotal_PLP\tin_v4.2.1\tin_v5.0q\tpct_of_gene_lost\n"
for g in PMS2 PKD1 NPHP1 ELANE OTOA KANSL1 FOXP1 RP1L1 CHAMP1 GJA8 \
         DLL3 MCM3AP TUBB8 GATA4 SDHA TNNI3 MYO7A MAPT SHANK3 RAC1 \
         ATP13A2 KCNQ2 HERC2 RNF216 VWA1; do
  awk -F'\t' -v g="$g" '$4 ~ "^"g":"' clinvar_plp.bed > .g.bed || true
  tot=$(wc -l < .g.bed)
  if [ "$tot" -gt 0 ]; then
    a=$(bedtools intersect -u -a .g.bed -b "$V421" | wc -l)
    b=$(bedtools intersect -u -a .g.bed -b "$V50Q" | wc -l)
    awk -v g="$g" -v t="$tot" -v a="$a" -v b="$b" \
      'BEGIN{ printf "%s\t%d\t%d\t%d\t%.1f%%\n", g, t, a, b, (a>0 ? 100*(a-b)/a : 0) }'
  fi
done
rm -f .g.bed

echo
echo "=== B. segmental duplication overlap ==="
wget -q -c $GIAB/genome-stratifications/v3.6/GRCh38@all/SegmentalDuplications/GRCh38_segdups.bed.gz
zcat GRCh38_segdups.bed.gz | grep -v '^#' \
  | awk 'BEGIN{OFS="\t"}{ c=$1; if (c !~ /^chr/) c="chr"c;
      if (c ~ /^chr([1-9]|1[0-9]|2[0-2])$/) print c,$2,$3 }' \
  | sort -k1,1 -k2,2n | bedtools merge -i - > segdups.bed

bpsum () { awk '{s+=$3-$2} END{printf "%d\n", s+0}' "$1"; }
lost_tot=$(bpsum lost_in_v50q.bed)
gain_tot=$(bpsum gained_in_v50q.bed)
lost_sd=$(bedtools intersect -a lost_in_v50q.bed -b segdups.bed | bpsum /dev/stdin)
gain_sd=$(bedtools intersect -a gained_in_v50q.bed -b segdups.bed | bpsum /dev/stdin)

awk -v lt=$lost_tot -v ls=$lost_sd -v gt=$gain_tot -v gs=$gain_sd 'BEGIN{
  printf "lost bases      : %d, of which %d (%.1f%%) in segdups\n", lt, ls, 100*ls/lt;
  printf "gained bases    : %d, of which %d (%.1f%%) in segdups\n", gt, gs, 100*gs/gt;
}'

echo
echo "=== C. ACMG SF v3.3 genes ==="
if [ -s acmg_sf_v33.txt ]; then
  bedtools intersect -u -a clinvar_plp.bed -b lost_in_v50q.bed \
    | cut -f4 | sed 's/|.*//;s/:.*//' | sort -u > .lostgenes.txt
  sort -u acmg_sf_v33.txt > .acmg.txt
  n_gene=$(comm -12 .lostgenes.txt .acmg.txt | wc -l)
  echo "distinct ACMG SF genes affected : $n_gene"
  comm -12 .lostgenes.txt .acmg.txt | tr '\n' ' '; echo
  rm -f .lostgenes.txt .acmg.txt
else
  echo "acmg_sf_v33.txt not found, skipping"
  echo "build it as one gene symbol per line from Lee K et al., Genet Med 2025;27(8):101454"
fi
