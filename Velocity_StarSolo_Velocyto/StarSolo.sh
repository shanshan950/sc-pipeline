#!/bin/bash
path=/mnt/rstor/genetics/JinLab/ssz20/zshanshan/differentiation_ATAC/processed/RNA_velocity/STARsolo/
# H1
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix H1 \
--readFilesIn SRR10902817_2.fastq SRR10902817_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/H1.barcode \
--soloFeatures Gene Velocyto
# S0
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S0 \
--readFilesIn SRR10902818_2.fastq SRR10902818_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S0.barcode \
--soloFeatures Gene Velocyto
# S1D1.barcode
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S1D1 \
--readFilesIn SRR10902819_2.fastq SRR10902819_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S1D1.barcode \
--soloFeatures Gene Velocyto
# S1D2.barcode
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S1D2 \
--readFilesIn SRR10902820_2.fastq SRR10902820_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S1D2.barcode \
--soloFeatures Gene Velocyto

# S2D1.barcode
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S2D1 \
--readFilesIn SRR10902821_2.fastq SRR10902821_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S2D1.barcode \
--soloFeatures Gene Velocyto

# S2D2.barcode
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S2D2 \
--readFilesIn SRR10902822_2.fastq SRR10902822_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S2D2.barcode \
--soloFeatures Gene Velocyto

# S2D3.barcode
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S2D3 \
--readFilesIn SRR10902823_2.fastq,SRR10902824_2.fastq SRR10902823_1.fastq.20base,SRR10902824_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S2D3.barcode \
--soloFeatures Gene Velocyto
# S3.barcode
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S3 \
--readFilesIn SRR10902825_2.fastq,SRR10902826_2.fastq SRR10902825_1.fastq.20base,SRR10902826_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S3.barcode \
--soloFeatures Gene Velocyto
# S4
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S4 \
--readFilesIn SRR10902827_2.fastq SRR10902827_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S4.barcode \
--soloFeatures Gene Velocyto
# S5
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S5 \
--readFilesIn SRR10902830_2.fastq SRR10902830_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S5.barcode \
--soloFeatures Gene Velocyto
# S6
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S6 \
--readFilesIn SRR10902832_2.fastq,SRR10902833_2.fastq SRR10902832_1.fastq.20base,SRR10902833_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S6.barcode \
--soloFeatures Gene Velocyto
# S7
$path/STAR \
--runThreadN 8 --genomeDir /mnt/rstor/genetics/JinLab/ssz20/reference/hg19_STARIndex/ --outFileNamePrefix S7 \
--readFilesIn SRR10902835_2.fastq SRR10902835_1.fastq.20base --soloBarcodeMate 0 --soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
--soloCBwhitelist $path/S7.barcode \
--soloFeatures Gene Velocyto

