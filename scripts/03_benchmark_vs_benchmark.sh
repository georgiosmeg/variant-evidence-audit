#!/usr/bin/env bash
set -euo pipefail

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
HG002=$GIAB/AshkenazimTrio/HG002_NA24385_son

V421BED=HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.clean.bed
V50QBED=HG002_GRCh38_v5.0q_smvar.benchmark.clean.bed
A_VCF=HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
B_VCF=HG002_GRCh38_v5.0q_smvar.vcf.gz

for f in "$V421BED" "$V50QBED"; do
  [ -s "$f" ] || { echo "missing $f, run 01 first" >&2; exit 1; }
done

echo "== downloading =="
wget -q -c $GIAB/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
wget -q -c $HG002/NISTv4.2.1/GRCh38/$A_VCF
wget -q -c $HG002/NISTv4.2.1/GRCh38/$A_VCF.tbi
wget -q -c $HG002/v5.0q/$B_VCF
wget -q -c $HG002/v5.0q/$B_VCF.tbi

REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
if ! samtools faidx "$REF" 2>/dev/null; then
  echo "== recompressing reference as bgzf =="
  zcat "$REF" | bgzip -@ 4 > GRCh38.fa.gz
  samtools faidx GRCh38.fa.gz
  REF=GRCh38.fa.gz
fi
echo "$REF" > .refpath

echo
echo "== diagnostics =="
for v in "$A_VCF" "$B_VCF"; do
  echo "--- $v"
  echo -n "  sample(s)    : "; bcftools query -l "$v" | tr '\n' ' '; echo
  echo -n "  total records: "; bcftools index -n "$v"
  echo -n "  first contig : "; bcftools view -h "$v" | grep -m1 '^##contig'
  echo    "  FILTER values in first 200k records:"
  bcftools view -H "$v" | head -200000 | cut -f7 | sort | uniq -c | sort -rn | head -5 | sed 's/^/    /'
done

echo
echo "== shared region =="
bedtools intersect -a "$V421BED" -b "$V50QBED" | sort -k1,1 -k2,2n | bedtools merge -i - > shared.bed
bpsum () { awk '{s+=$3-$2} END{printf "%d\n", s+0}' "$1"; }
echo "shared confident bases : $(bpsum shared.bed)"
echo "shared intervals       : $(wc -l < shared.bed)"

bcftools view -T shared.bed "$A_VCF" -Oz -o A.raw.vcf.gz
bcftools view -T shared.bed "$B_VCF" -Oz -o B.raw.vcf.gz
bcftools index -f -t A.raw.vcf.gz
bcftools index -f -t B.raw.vcf.gz

nA=$(bcftools index -n A.raw.vcf.gz)
nB=$(bcftools index -n B.raw.vcf.gz)
echo "v4.2.1 records in shared region : $nA"
echo "v5.0q  records in shared region : $nB"
[ "$nB" -gt 0 ] || { echo "v5.0q returned nothing, check contig naming and the index" >&2; exit 1; }

echo
echo "== L1: raw, no normalisation =="
rm -rf isec_L1; bcftools isec -p isec_L1 A.raw.vcf.gz B.raw.vcf.gz
L1a=$(grep -vc '^#' isec_L1/0000.vcf || true)
L1b=$(grep -vc '^#' isec_L1/0001.vcf || true)
L1s=$(grep -vc '^#' isec_L1/0002.vcf || true)
awk -v a=$L1a -v b=$L1b -v s=$L1s 'BEGIN{
  printf "  v4.2.1 only %d, v5.0q only %d, shared %d\n", a, b, s;
  printf "  disagreement rate: %.2f%%\n", 100*(a+b)/(a+b+s) }'

echo
echo "== L2: left-aligned, atomised, duplicates dropped =="
ATOM=""
bcftools norm --help 2>&1 | grep -q -- '--atomize' && ATOM="-a"
for s in A B; do
  bcftools norm -f "$REF" $ATOM -m -any -c w "$s.raw.vcf.gz" 2>/dev/null \
    | bcftools norm -d exact -Oz -o "$s.norm.vcf.gz"
  bcftools index -f -t "$s.norm.vcf.gz"
done
rm -rf isec_L2; bcftools isec -p isec_L2 A.norm.vcf.gz B.norm.vcf.gz
L2a=$(grep -vc '^#' isec_L2/0000.vcf || true)
L2b=$(grep -vc '^#' isec_L2/0001.vcf || true)
L2s=$(grep -vc '^#' isec_L2/0002.vcf || true)
awk -v a=$L2a -v b=$L2b -v s=$L2s 'BEGIN{
  printf "  v4.2.1 only %d, v5.0q only %d, shared %d\n", a, b, s;
  printf "  disagreement rate: %.2f%%\n", 100*(a+b)/(a+b+s) }'

cat isec_L2/0000.vcf isec_L2/0001.vcf | grep -v '^#' \
  | awk 'BEGIN{OFS="\t"}{ print $1, $2-1, $2-1+length($4) }' \
  | sort -k1,1 -k2,2n -u > disagreements.bed
grep -v '^#' isec_L2/0002.vcf \
  | awk 'BEGIN{OFS="\t"}{ print $1, $2-1, $2-1+length($4) }' \
  | sort -k1,1 -k2,2n -u > concordant.bed

echo
echo "distinct disagreeing positions : $(wc -l < disagreements.bed)"
echo "distinct concordant positions  : $(wc -l < concordant.bed)"
