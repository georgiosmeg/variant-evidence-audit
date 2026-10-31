#!/usr/bin/env bash
set -uo pipefail

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
V50QBED=HG002_GRCh38_v5.0q_smvar.benchmark.clean.bed
B_VCF=HG002_GRCh38_v5.0q_smvar.vcf.gz
REF=$(cat .refpath 2>/dev/null || echo GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz)

for f in isec_L2/0000.vcf isec_L2/0001.vcf isec_L2/0002.vcf \
         clinvar_plp.bed shared.bed segdups.bed disagreements.bed concordant.bed \
         A.norm.vcf.gz B.norm.vcf.gz B.raw.vcf.gz "$REF" "$B_VCF"; do
  [ -s "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

cntgz () { local n; n=$(zgrep -vc '^#' "$1" 2>/dev/null || true); echo "${n:-0}"; }

echo "=== 1. denominator for the clinical overlap ==="
cat concordant.bed disagreements.bed | sort -k1,1 -k2,2n -u > allbench.bed
clin_total=$(wc -l < clinvar_plp.bed)
clin_on_bench=$(bedtools intersect -u -a clinvar_plp.bed -b allbench.bed | wc -l)
clin_dis=$(bedtools intersect -u -a clinvar_plp.bed -b disagreements.bed | wc -l)
echo "ClinVar P/LP positions, autosomes               : $clin_total"
echo "...that coincide with a benchmark record at all : $clin_on_bench"
echo "...where the two benchmarks disagree            : $clin_dis"
awk -v a=$clin_dis -v b=$clin_on_bench 'BEGIN{
  if (b>0) printf "disagreement rate among those       : %.2f%% (%d of %d)\n", 100*a/b, a, b }'

echo
echo "=== 2. variant class of the disagreements ==="
for side in "v4.2.1-only:isec_L2/0000.vcf" "v5.0q-only:isec_L2/0001.vcf"; do
  label=${side%%:*}; file=${side##*:}
  echo "--- $label"
  grep -v '^#' "$file" \
    | awk '{ if (length($4)==1 && length($5)==1) c="SNV";
             else if (length($4) > length($5)) c="DEL";
             else if (length($4) < length($5)) c="INS";
             else c="MNV"; print c }' \
    | sort | uniq -c | sort -rn | sed 's/^/    /'
done

echo
echo "=== 3. FILTER audit ==="
echo "--- inside the shared benchmark region"
bcftools view -H B.raw.vcf.gz | cut -f7 | sort | uniq -c | sort -rn | head -8 | sed 's/^/    /'
echo "--- whole v5.0q VCF, genome-wide"
bcftools view -H "$B_VCF" | cut -f7 | sort | uniq -c | sort -rn | head -8 | sed 's/^/    /'
echo "--- inside the v5.0q benchmark BED only"
bcftools view -H -T "$V50QBED" "$B_VCF" | cut -f7 | sort | uniq -c | sort -rn | head -8 | sed 's/^/    /'

echo
echo "=== 4. PMS2 locus check ==="
awk -F'\t' '$4 ~ /^PMS2:/' clinvar_plp.bed | sort -k1,1 -k2,2n | bedtools merge -i - > pms2_span.bed
awk '{s+=$3-$2} END{printf "PMS2 span: %d bp across %d intervals\n", s, NR}' pms2_span.bed
ov=$(bedtools intersect -a pms2_span.bed -b "$V50QBED" | awk '{s+=$3-$2} END{print s+0}')
echo "overlap with the v5.0q benchmark region: $ov bp"

echo
echo "=== 5. stratification files ==="
wget -q -c $GIAB/genome-stratifications/v3.6/GRCh38@all/Union/GRCh38_alllowmapandsegdupregions.bed.gz
wget -q -c $GIAB/genome-stratifications/v3.6/GRCh38@all/LowComplexity/GRCh38_AllTandemRepeats.bed.gz
prep () {
  zcat "$1" | grep -v '^#' \
    | awk 'BEGIN{OFS="\t"}{ c=$1; if (c !~ /^chr/) c="chr"c;
        if (c ~ /^chr([1-9]|1[0-9]|2[0-2])$/) print c,$2,$3 }' \
    | sort -k1,1 -k2,2n | bedtools merge -i - > "$2"
}
[ -s hard.bed ] || prep GRCh38_alllowmapandsegdupregions.bed.gz hard.bed
[ -s tr.bed ]   || prep GRCh38_AllTandemRepeats.bed.gz tr.bed
echo "hard.bed and tr.bed ready"

echo
echo "=== 6. L3 haplotype-aware ==="
jmaj () {
  local s
  s=$("$1" -version 2>&1 | head -1 | sed -E 's/^[^"]*"([^"]+)".*$/\1/') || return 1
  case "$s" in 1.*) echo "$s" | cut -d. -f2 ;; *) echo "$s" | cut -d. -f1 ;; esac
}
JAVA=""
for cand in "$(command -v java 2>/dev/null)" "${CONDA_PREFIX:-/nonexistent}/bin/java" \
            /usr/lib/jvm/*/bin/java /opt/*/bin/java; do
  [ -n "$cand" ] && [ -x "$cand" ] || continue
  m=$(jmaj "$cand" 2>/dev/null)
  [ -n "$m" ] && [ "$m" -ge 11 ] 2>/dev/null && { JAVA="$cand"; break; }
done
if [ -z "$JAVA" ]; then
  echo 'no Java 11 or newer found, run: conda install -c conda-forge "openjdk>=11"'
  exit 0
fi
echo "java: $JAVA ($(jmaj "$JAVA"))"

RTGJAR=$(find "${CONDA_PREFIX:-/nonexistent}" -name RTG.jar 2>/dev/null | head -1)
[ -n "$RTGJAR" ] || { echo "RTG.jar not found" >&2; exit 1; }
echo "jar : $RTGJAR"

if [ ! -d GRCh38.sdf ]; then
  "$JAVA" -Xmx12g -jar "$RTGJAR" format -o GRCh38.sdf "$REF" \
    || { zcat "$REF" > GRCh38.fa; "$JAVA" -Xmx12g -jar "$RTGJAR" format -o GRCh38.sdf GRCh38.fa; }
fi

rm -rf vcfeval_out
"$JAVA" -Xmx12g -jar "$RTGJAR" vcfeval \
  -b A.norm.vcf.gz -c B.norm.vcf.gz -e shared.bed \
  -t GRCh38.sdf -o vcfeval_out --ref-overlap --all-records

L3a=$(cntgz vcfeval_out/fn.vcf.gz)
L3b=$(cntgz vcfeval_out/fp.vcf.gz)
L3s=$(cntgz vcfeval_out/tp.vcf.gz)
echo
echo "v4.2.1 only : $L3a"
echo "v5.0q  only : $L3b"
echo "concordant  : $L3s"
awk -v a=$L3a -v b=$L3b -v s=$L3s 'BEGIN{
  if (a+b+s>0) printf "disagreement rate: %.3f%%\n", 100*(a+b)/(a+b+s) }'
L1=$(( $(grep -vc '^#' isec_L1/0000.vcf || true) + $(grep -vc '^#' isec_L1/0001.vcf || true) ))
awk -v l1=$L1 -v l3=$((L3a+L3b)) 'BEGIN{
  if (l1>0) printf "removed L1 to L3 : %.1f%%\n", 100*(l1-l3)/l1 }'

echo
echo "=== 7. surviving disagreements ==="
zcat vcfeval_out/fn.vcf.gz vcfeval_out/fp.vcf.gz 2>/dev/null | grep -v '^#' \
  | awk 'BEGIN{OFS="\t"}{ print $1, $2-1, $2-1+length($4) }' \
  | sort -k1,1 -k2,2n -u > real_disagreements.bed
n_real=$(wc -l < real_disagreements.bed)
n_con=$(wc -l < concordant.bed)
echo "distinct surviving positions: $n_real"

printf "\n%-10s %12s %12s %10s\n" "stratum" "surviving" "concordant" "enrichment"
for s in segdups hard tr; do
  d=$(bedtools intersect -u -a real_disagreements.bed -b $s.bed | wc -l)
  c=$(bedtools intersect -u -a concordant.bed -b $s.bed | wc -l)
  awk -v n="$s" -v d=$d -v c=$c -v nd=$n_real -v nc=$n_con 'BEGIN{
    if (nd>0 && nc>0) { pd=100*d/nd; pc=100*c/nc;
      printf "%-10s %11.1f%% %11.1f%% %9.1fx\n", n, pd, pc, (pc>0?pd/pc:0) } }'
done

echo
echo "ClinVar P/LP positions among the survivors:"
bedtools intersect -u -a clinvar_plp.bed -b real_disagreements.bed | wc -l
bedtools intersect -u -a clinvar_plp.bed -b real_disagreements.bed \
  | cut -f4 | sed 's/|.*//' | sort | uniq -c | sort -rn | head -15 \
  | tee genes_real_disagreement.txt
