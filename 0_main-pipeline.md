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