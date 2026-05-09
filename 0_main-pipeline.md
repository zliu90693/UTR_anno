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
mkdir -p fastqc_out
fastqc unmapped.fastq -o fastqc_out/
```