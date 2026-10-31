#!/usr/bin/env bash
set -euo pipefail

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son
CV=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38

for tool in bedtools bcftools wget awk sort; do
  command -v "$tool" >/dev/null 2>&1 || { echo "missing: $tool" >&2; exit 1; }
done

echo "== downloading =="
wget -q -c $GIAB/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
wget -q -c $GIAB/v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed
wget -q -c $GIAB/v5.0q/HG002_CHM13v2.0_v5.0q_smvar.benchmark.bed
wget -q -c $CV/clinvar.vcf.gz
wget -q -c $CV/clinvar.vcf.gz.tbi

echo "== normalising regions =="
for f in HG002_*.bed; do
  case "$f" in *.clean.bed) continue;; esac
  grep -v '^#' "$f" \
    | awk 'BEGIN{OFS="\t"}{ c=$1; if (c !~ /^chr/) c="chr"c;
        if (c ~ /^chr([1-9]|1[0-9]|2[0-2])$/) print c,$2,$3 }' \
    | sort -k1,1 -k2,2n | bedtools merge -i - > "${f%.bed}.clean.bed"
done

V421=HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.clean.bed
V50Q=HG002_GRCh38_v5.0q_smvar.benchmark.clean.bed
V50C=HG002_CHM13v2.0_v5.0q_smvar.benchmark.clean.bed

bp () { awk '{s+=$3-$2} END{printf "%d\n", s+0}' "$1"; }

echo
echo "=== evaluable bases, autosomes ==="
echo "v4.2.1 on GRCh38 : $(bp $V421)"
echo "v5.0q  on GRCh38 : $(bp $V50Q)"
echo "v5.0q  on CHM13  : $(bp $V50C)"

bedtools subtract -a $V50Q -b $V421 > gained_in_v50q.bed
bedtools subtract -a $V421 -b $V50Q > lost_in_v50q.bed

echo
echo "bases gained in v5.0q : $(bp gained_in_v50q.bed)"
echo "bases lost   in v5.0q : $(bp lost_in_v50q.bed)"

echo
echo "=== per chromosome ==="
printf "chrom\tv4.2.1_bp\tv5.0q_bp\tdelta\n"
for c in $(seq 1 22); do
  a=$(awk -v c=chr$c '$1==c{s+=$3-$2} END{print s+0}' $V421)
  b=$(awk -v c=chr$c '$1==c{s+=$3-$2} END{print s+0}' $V50Q)
  printf "chr%s\t%s\t%s\t%s\n" "$c" "$a" "$b" "$((b-a))"
done

echo
echo "== extracting ClinVar pathogenic / likely pathogenic =="
bcftools view -i 'INFO/CLNSIG ~ "Pathogenic" || INFO/CLNSIG ~ "Likely_pathogenic"' \
  clinvar.vcf.gz -Oz -o clinvar_plp.vcf.gz
bcftools index -f -t clinvar_plp.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%INFO/GENEINFO\t%INFO/CLNSIG\n' clinvar_plp.vcf.gz \
  | awk 'BEGIN{OFS="\t"}{ c=$1; if (c !~ /^chr/) c="chr"c;
      if (c ~ /^chr([1-9]|1[0-9]|2[0-2])$/) print c, $2-1, $2-1+length($3), $4, $5 }' \
  | sort -k1,1 -k2,2n > clinvar_plp.bed

total=$(wc -l < clinvar_plp.bed)
in421=$(bedtools intersect -u -a clinvar_plp.bed -b $V421 | wc -l)
in50q=$(bedtools intersect -u -a clinvar_plp.bed -b $V50Q | wc -l)
lost=$(bedtools intersect -u -a clinvar_plp.bed -b lost_in_v50q.bed | wc -l)
gained=$(bedtools intersect -u -a clinvar_plp.bed -b gained_in_v50q.bed | wc -l)

echo
echo "=== ClinVar pathogenic / likely pathogenic, autosomes ==="
echo "total                               : $total"
echo "inside v4.2.1 evaluable region      : $in421"
echo "inside v5.0q  evaluable region      : $in50q"
echo "became UNEVALUABLE between releases : $lost"
echo "became EVALUABLE between releases   : $gained"

bedtools intersect -u -a clinvar_plp.bed -b lost_in_v50q.bed \
  | cut -f4 | sed 's/|.*//' | sort | uniq -c | sort -rn | head -25 \
  > genes_losing_coverage.txt

echo
echo "=== top genes losing benchmark coverage ==="
cat genes_losing_coverage.txt
