```bash
conda activate fastq2matrix
./3_cellranger-count.sh -p Jones_NEE_2023_Lzep -l 20 -t L_zep_mkref_out -c true
```

```bash
conda activate UTR_anno
samtools merge \
    -@ 8 \
    merged.bam \
    "./Jones_NEE_2023_Lzep/cellranger-count-out/LZEP-Queen/outs/possorted_genome_bam.bam" \
    "./Jones_NEE_2023_Lzep/cellranger-count-out/LZEP-Worker/outs/possorted_genome_bam.bam"

samtools index merged.bam
```

```bash
mv merged.bam "./Jones_NEE_2023_Lzep/ref"
mv merged.bam.bai "./Jones_NEE_2023_Lzep/ref"
```

```bash
peaks2utr "./Jones_NEE_2023_Lzep/ref/Lasioglossum_zephyrus.gtf" merged.bam -p 15 --gtf -o Lasioglossum_zephyrus_queen_fixed.gtf
```
[check 3'utr](./ref-inspection.ipynb)
```bash
# mv Lasioglossum_zephyrus_queen_fixed.gtf ./Jones_NEE_2023_Lzep/ref/
./2_make-ref.sh -p Jones_NEE_2023_Lzep -g Lasioglossum_zephyrus_queen_fixed.gtf -f Lasioglossum_zephyrus_queen_fixed.filtered.gtf -r Lasioglossum_zephyrus.fasta -m L_zep_mkref_fixed_out 
```
```bash
./3_cellranger-count.sh -p Jones_NEE_2023_Lzep -l 20 -t L_zep_mkref_fixed_out -c false
```
```bash
samtools flagstat ./Jones_NEE_2023_Lzep/cellranger-count-out/LZEP-Queen-old/outs/possorted_genome_bam.bam

samtools view -f 4 Jones_NEE_2023_Lzep/cellranger-count-out/LZEP-Queen-old/outs/possorted_genome_bam.bam | head -1000 | \
awk '{print ">"$1"\n"$10}' > unmapped_sample.fa

samtools idxstats Jones_NEE_2023_Lzep/cellranger-count-out/LZEP-Queen-old/outs/possorted_genome_bam.bam | sort -k3 -rn | head -100 # 发现三个unplaced的reads极高, blast发现是rRNA, 且三个基本在GTF中无注释

# 确定了已map到基因组中的reads中, 有相当一部分map到了rRNA上

samtools faidx Jones_NEE_2023_Lzep/ref/Lasioglossum_zephyrus.fasta LZEP_unplaced_8928 > ./metadata/unplaced_8928.fa
```
但是, 为什么有80%的reads没有map到基因组上, 原因仍不确定:
```bash
samtools fastq -f 4 Jones_NEE_2023_Lzep/cellranger-count-out/LZEP-Queen-old/outs/possorted_genome_bam.bam > unmapped.fastq
```
```bash
mkdir -p fastqc_out/fastqc_out_raw
fastqc unmapped.fastq -o fastqc_out/fastqc_out_raw
```

```bash
echo "Queen R1 reads:"
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Queen/LZEP-Queen_S1_L001_R1_001.fastq.gz | wc -l | awk '{print $1/4}'
echo "Queen R2 reads:"
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Queen/LZEP-Queen_S1_L001_R2_001.fastq.gz  | wc -l | awk '{print $1/4}'
```
Queen R1 reads: 264430935, Queen R2 reads: 264430935
```bash
fastp \
    -i Jones_NEE_2023_Lzep/fastq/LZEP-Queen/LZEP-Queen_S1_L001_R1_001.fastq.gz \
    -I Jones_NEE_2023_Lzep/fastq/LZEP-Queen/LZEP-Queen_S1_L001_R2_001.fastq.gz \
    -o Jones_NEE_2023_Lzep/fastq/LZEP-Queen/trimed_R1.fastq.gz \
    -O Jones_NEE_2023_Lzep/fastq/LZEP-Queen/trimed_R2.fastq.gz \
    --thread 16
```
```text
Read1 before filtering:
total reads: 264430935
total bases: 7404066180
Q20 bases: 7403879504(99.9975%)
Q30 bases: 7403879504(99.9975%)

Read1 after filtering:
total reads: 262982755
total bases: 7363517140
Q20 bases: 7363517140(100%)
Q30 bases: 7363517140(100%)

Read2 before filtering:
total reads: 264430935
total bases: 24856507890
Q20 bases: 24721937208(99.4586%)
Q30 bases: 24721937208(99.4586%)

Read2 aftering filtering:
total reads: 262982755
total bases: 21828365098
Q20 bases: 21828365098(100%)
Q30 bases: 21828365098(100%)

Filtering result:
reads passed filter: 525965510
reads failed due to low quality: 2875354
reads failed due to too many N: 21006
reads failed due to too short: 0
reads with adapter trimmed: 87675300
bases trimmed due to adapters: 2893284900
```

```bash
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Queen/trimed_R1.fastq.gz | wc -l | awk '{print $1/4}'
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Queen/trimed_R2.fastq.gz  | wc -l | awk '{print $1/4}'
```
262982755, 262982755
```bash
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Worker/LZEP-Worker_S1_L001_R1_001.fastq.gz | wc -l | awk '{print $1/4}'
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Worker/LZEP-Worker_S1_L001_R2_001.fastq.gz  | wc -l | awk '{print $1/4}'
# 169559757
# 169559757
fastp \
    -i Jones_NEE_2023_Lzep/fastq/LZEP-Worker/LZEP-Worker_S1_L001_R1_001.fastq.gz \
    -I Jones_NEE_2023_Lzep/fastq/LZEP-Worker/LZEP-Worker_S1_L001_R2_001.fastq.gz \
    -o Jones_NEE_2023_Lzep/fastq/LZEP-Worker/trimed_R1.fastq.gz \
    -O Jones_NEE_2023_Lzep/fastq/LZEP-Worker/trimed_R2.fastq.gz \
    --thread 16

# Read1 before filtering:
# total reads: 169559757
# total bases: 4747673196
# Q20 bases: 4747564836(99.9977%)
# Q30 bases: 4747564836(99.9977%)
# Q40 bases: 0(0%)

# Read2 before filtering:
# total reads: 169559757
# total bases: 15938617158
# Q20 bases: 15886433904(99.6726%)
# Q30 bases: 15886433904(99.6726%)
# Q40 bases: 0(0%)

# Read1 after filtering:
# total reads: 168998724
# total bases: 4731964272
# Q20 bases: 4731964272(100%)
# Q30 bases: 4731964272(100%)
# Q40 bases: 0(0%)

# Read2 after filtering:
# total reads: 168998724
# total bases: 13821275316
# Q20 bases: 13821275316(100%)
# Q30 bases: 13821275316(100%)
# Q40 bases: 0(0%)

# Filtering result:
# reads passed filter: 337997448
# reads failed due to low quality: 1117324
# reads failed due to too many N: 4742
# reads failed due to too short: 0
# reads failed due to adapter dimer: 0
# reads with adapter trimmed: 62618588
# bases trimmed due to adapters: 2066413404

# Duplication rate: 19.267%

# Insert size peak (evaluated by paired-end reads): 28

# JSON report: fastp.json
# HTML report: fastp.html

zcat Jones_NEE_2023_Lzep/fastq/LZEP-Worker/LZEP-Worker_S1_L001_R1_001.fastq.gz | wc -l | awk '{print $1/4}'
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Worker/LZEP-Worker_S1_L001_R2_001.fastq.gz  | wc -l | awk '{print $1/4}'
# 169559757
# 169559757
```

```bash
mkdir -p fastqc_out/fastqc_out_filtered_queen
fastqc Jones_NEE_2023_Lzep/fastq/LZEP-Queen/trimed_R2.fastq.gz -o fastqc_out/fastqc_out_filtered_queen
```
```bash
mkdir -p fastqc_out/fastqc_out_filtered_worker
fastqc Jones_NEE_2023_Lzep/fastq/LZEP-Worker/trimed_R2.fastq.gz -o fastqc_out/fastqc_out_filtered_worker
```
检查fastp是否对R1造成了破坏:
```bash
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Queen/trimed_R1.fastq.gz | awk 'NR%4==2 {print length}' | sort -n | uniq -c | head -20
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Queen/trimed_R2.fastq.gz | awk 'NR%4==2 {print length}' | sort -n | uniq -c | head -20

zcat Jones_NEE_2023_Lzep/fastq/LZEP-Worker/trimed_R1.fastq.gz | awk 'NR%4==2 {print length}' | sort -n | uniq -c | head -20
zcat Jones_NEE_2023_Lzep/fastq/LZEP-Worker/trimed_R2.fastq.gz | awk 'NR%4==2 {print length}' | sort -n | uniq -c | head -20
```