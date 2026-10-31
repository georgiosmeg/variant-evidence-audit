# variant-evidence-audit

What happens to the measurement when a genomics benchmark changes?

Short answer: the two most recent GIAB benchmarks for the same sample on the same reference genome agree almost perfectly about everything they both look at. The problem is what the newer one stopped looking at. 2,670 clinically annotated pathogenic positions left the measurable region, and thirteen medically important genes went to zero coverage, including all of *PMS2*.

This repository holds the scripts, the logs and the small derived files behind that result. Everything runs on public data.

---

## Why this matters

AI models do most of the variant calling in clinical genomics. Their reported accuracy is measured against a benchmark, the benchmark is tied to one reference genome, and the reference genome is changing. If a benchmark update quietly removes the hardest clinically relevant cases from the test set, every model scored on it looks better and no aggregate accuracy number says why.

That is the failure mode this analysis documents. It is not a criticism of GIAB, whose exclusions are deliberate and well reasoned. It is a gap downstream: nothing propagates those exclusions into the numbers people compare.

---

## Findings

### 1. What left the test set

Comparing the HG002 v4.2.1 and v5.0q benchmark regions on GRCh38, autosomes only.

| | Count |
| --- | ---: |
| ClinVar pathogenic / likely pathogenic positions | 339,610 |
| Measurable under v4.2.1 | 332,426 |
| Measurable under v5.0q | 336,011 |
| Became measurable | 6,930 |
| **Stopped being measurable** | **2,670** |
| Share of the loss in the top 25 genes | 87% |
| Share of the loss in *PMS2* alone | 33% |

The net change is positive. Where the losses land is the point.

### 2. Thirteen genes go to zero

| Gene | P/LP positions | In v4.2.1 | In v5.0q |
| --- | ---: | ---: | ---: |
| PMS2 | 935 | 885 | 0 |
| NPHP1 | 146 | 146 | 0 |
| ELANE | 119 | 119 | 0 |
| OTOA | 95 | 95 | 0 |
| KANSL1 | 132 | 91 | 0 |
| RP1L1 | 82 | 77 | 0 |
| CHAMP1 | 61 | 61 | 0 |
| GJA8 | 60 | 60 | 0 |
| TUBB8 | 46 | 42 | 0 |
| GATA4 | 41 | 41 | 0 |
| MAPT | 36 | 36 | 0 |
| RAC1 | 33 | 33 | 0 |
| VWA1 | 22 | 22 | 0 |

Full list including partial losses in `results/genes_losing_coverage.txt`.

*PMS2* causes Lynch syndrome. It is hard to call because its pseudogene *PMS2CL* carries near-identical sequence across several exons inside a segmental duplication, and bringing it into the benchmark was a flagship result of v4.2.1, which took coverage from 25.9% in v3.3.2 to 85.6%. The exclusion in v5.0q is the whole locus, not just those positions: the 35,715 bp span containing *PMS2*'s clinically annotated variants has zero overlap with the v5.0q benchmark region.

Note that *PKD1* appears in the top-25 loss list but is net positive overall, 2,000 positions covered under v4.2.1 against 2,337 under v5.0q. It lost some territory and gained more. Describe it as mixed.

### 3. The territory that was traded

| | Bases | In segmental duplications |
| --- | ---: | ---: |
| Given up by v5.0q | 46,035,747 | 28,557,338 (62.0%) |
| Gained by v5.0q | 83,105,444 | 11,741,212 (14.1%) |

A 4.4-fold enrichment. Easy new ground for hard old ground.

One number to be careful with: the v5.0q preprint reports 9.1% more reference coverage. On the autosomes, where like-for-like comparison is possible, the gain is 1.5%, from 2,542,242,843 bases to 2,579,312,540. Most of the headline figure is the sex chromosomes, which v4.2.1 never covered at all.

### 4. Where both benchmarks look, do they agree?

Shared high-confidence territory is 2,496,207,096 bases across 445,112 intervals. v4.2.1 has 3,811,741 records in it, v5.0q has 3,813,687.

| Comparison | v4.2.1 only | v5.0q only | Total | Rate |
| --- | ---: | ---: | ---: | ---: |
| L1, raw records | 5,762 | 7,708 | 13,470 | 0.35% |
| L2, left-aligned and atomised | 5,534 | 7,611 | 13,145 | 0.34% |
| **L3, haplotype-aware** | **145** | **1,152** | **1,297** | **0.034%** |

**90.4% of the apparent disagreement was bookkeeping.** The split between rungs is the useful part. Normalisation, which is what most people actually run, removed 2.4%. Haplotype-aware comparison removed the other 88%.

### 5. What survives

| Stratum | Surviving disagreements | Concordant positions | Enrichment |
| --- | ---: | ---: | ---: |
| Segmental duplications | 83.6% | 2.4% | 35.4x |
| Low mappability or segdup | 86.6% | 6.1% | 14.1x |
| Tandem repeats | 8.2% | 4.9% | 1.7x |

The tandem-repeat signal was 43% of disagreements at L2 and collapses to 8.2% at L3, because those were the representation artefacts. What is left is segmental duplication, and it is lopsided: 1,152 positions claimed by the assembly-based benchmark against 145 by the mapping-based one. The 1,209 distinct surviving positions are in `results/real_disagreements.bed`.

### 6. No clinical disagreement, only clinical absence

Only 1,704 of the 339,610 ClinVar pathogenic positions coincide with a benchmark record at all, because HG002 is one person who does not carry most pathogenic variants. Of those 1,704, the benchmarks disagree about 4 at the naive level and none after haplotype-aware comparison.

So this is not two truth sets contradicting each other about patient-relevant variants. It is one of them no longer asking.

### 7. The benchmark flags its own uncertainty, then excludes it

| Scope | Unflagged | Flagged |
| --- | ---: | ---: |
| v5.0q genome-wide | 5,201,459 | 744,066 |
| Inside the v5.0q benchmark region | 4,404,375 | 0 |
| Inside the shared region | 3,813,687 | 0 |

744,066 records carry a `HET1`, `HET2`, `GAP1` or `GAP2` flag marking something the assembly could not resolve, and not one sits inside the benchmark region. GIAB flags what it could not resolve and then excludes it. That is the mechanism behind finding 2.

---

## Reproducing

Dependencies: `bcftools`, `bedtools`, `samtools`, `rtg-tools`, and a **Java 11 or newer** runtime.

```bash
conda install -c bioconda bcftools bedtools samtools rtg-tools
conda install -c conda-forge "openjdk>=11"

mkdir prelim && cd prelim
bash ../scripts/01_benchmark_churn_prelim.sh    # findings 1 and 2
bash ../scripts/02_followup_analyses.sh         # findings 2 and 3, per-gene detail
bash ../scripts/03_benchmark_vs_benchmark.sh    # finding 4, rungs L1 and L2
bash ../scripts/04_run_now.sh                   # finding 4 rung L3, findings 5 to 7
```

About 4 GB of downloads and a few hours end to end. No specialised compute.

### Two gotchas worth knowing

**rtg-tools 3.13 needs Java 11.** Conda installs OpenJDK 8 alongside it by default, which fails with `UnsupportedClassVersionError: class file version 55.0` against a runtime that reads up to 52.0. Install a newer JDK.

**On aarch64, the `rtg` launcher refuses to run** with `Sorry, you must be running a 64bit operating system`. That is the wrapper script checking `uname -m` against a list that predates ARM. The JVM is fine. Call the jar directly:

```bash
RTGJAR=$(find "$CONDA_PREFIX" -name RTG.jar | head -1)
java -Xmx12g -jar "$RTGJAR" vcfeval ...
```

---

## Caveats

- **v5.0q is a draft benchmark** and its exclusions are almost certainly deliberate. Nothing here says GIAB made a mistake. What is new is the clinical cost, in named genes and counted positions.
- **This is benchmark coverage and content, not model accuracy.** No caller was run. The analysis measures where truth exists and what it says, which is upstream of accuracy.
- **ClinVar positions are not HG002 genotypes.** The clinical denominator is 1,704, not 339,610.
- **One region was dropped.** `rtg vcfeval` abandoned `chr6:18111262-18111288` after 232,204 unresolved iterations and excluded it from results.
- **Autosomes only**, because v4.2.1 covers chromosomes 1 to 22 and nothing else.
- **No quality scores**, so vcfeval produced un-thresholded statistics only. This does not affect any count above.

---

## Data sources

All public, none redistributed here beyond small derived interval files.

| Source | Used for |
| --- | --- |
| [GIAB HG002 benchmarks](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/) | v4.2.1 and v5.0q regions and variant calls |
| [GIAB genome stratifications v3.6](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/) | Segmental duplication, low mappability and tandem repeat regions |
| [ClinVar GRCh38 VCF](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) | Pathogenic and likely pathogenic positions |
| [GRCh38 no-alt analysis set](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/) | Reference for normalisation and vcfeval |

Derived interval files in `results/` are computed from the above. Those sources carry their own terms; check them before redistributing anything larger.

---

## Contact

Georgios Megalovasilis, georgios.me@gmail.com

Work supported in part by a BlueDot Impact Rapid Grant application. See `docs/` for the proposal summary.
