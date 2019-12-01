#!/bin/bash
usage="./dropseq_pipline.sh <name prefix> <genome> <Reads1.fq> <Read2.fq> <estimated cell number>\n This Program run from fq raw data toward DGE"
########################
name=$1
genome=$2
fq1=$3
fq2=$4
cellnumber=20000
cellnumber2=40000
cellnumber5=100000

mkdir process_data/$name
myname=process_data/$name/$name
# requirement picard STAR
picard=software/picard-tools-1.93
star=software/STAR-2.5.1b
dropseq=software/Drop-seq_tools-1.0
genomeDir=reference/$genome
mylib=lib
##########################
$mylib/fastqmerge.pl $fq1 $fq2 20 50 > $myname.R1R2.fastq
java -jar -Xmx8g  $picard/FastqToSam.jar  F1=$myname.R1R2.fastq  O=$myname.R1R2.sam  SM=$myname
samtools view -h -o $myname.bam  $myname.R1R2.sam
$dropseq/TagBamWithReadSequenceExtended INPUT=$myname.bam OUTPUT=$myname.unaligned_tagged_Cell.bam SUMMARY=$myname.unaligned_tagged_Cell.bam_summary.txt BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1
$dropseq/TagBamWithReadSequenceExtended INPUT=$myname.unaligned_tagged_Cell.bam  OUTPUT=$myname.unaligned_tagged_CellMolecular.bam SUMMARY=$myname.unaligned_tagged_Molecular.bam_summary.txt BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1
$dropseq/FilterBAM TAG_REJECT=XQ INPUT=$myname.unaligned_tagged_CellMolecular.bam OUTPUT=$myname.unaligned_tagged_filtered.bam
$dropseq/TrimStartingSequence  INPUT=$myname.unaligned_tagged_filtered.bam  OUTPUT=$myname.unaligned_tagged_trimmed_smart.bam OUTPUT_SUMMARY=$myname.adapter_trimming_report.txt SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5
$dropseq/PolyATrimmer INPUT=$myname.unaligned_tagged_trimmed_smart.bam OUTPUT=$myname.unaligned_mc_tagged_polyA_filtered.bam OUTPUT_SUMMARY=$myname.polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6
samtools view $myname.unaligned_mc_tagged_polyA_filtered.bam |$mylib/trimsam.pl | samtools view -h -o $myname.unaligned_mc_tagged_polyA_filtered.trimed.bam -
samtools view -H $myname.unaligned_mc_tagged_polyA_filtered.bam > $myname.unaligned.header
samtools view -h $myname.unaligned_mc_tagged_polyA_filtered.trimed.bam | cat $myname.unaligned.header - > $myname.unaligned_mc_tagged_polyA_filtered.trimed.h.bam
java -jar -Xmx8g  $picard/SamToFastq.jar  I=$myname.unaligned_mc_tagged_polyA_filtered.trimed.h.bam F=$myname.unaligned_mc_tagged_polyA_filtered.trimed.fastq
$star/STAR --runThreadN 8 --genomeDir $genomeDir/${genome}STARindex --readFilesIn $myname.unaligned_mc_tagged_polyA_filtered.trimed.fastq --outFileNamePrefix $myname.star
java -Xmx8g -jar $picard/SortSam.jar I=$myname.starAligned.out.sam  O=$myname.aligned.sorted.bam SO=queryname
java -Xmx8g -jar $picard/MergeBamAlignment.jar UNMAPPED_BAM=$myname.unaligned_mc_tagged_polyA_filtered.trimed.h.bam ALIGNED_BAM=$myname.aligned.sorted.bam OUTPUT=$myname.merged.bam  REFERENCE_SEQUENCE=$genomeDir/$genome.fa INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false
$dropseq/TagReadWithGeneExon I=$myname.merged.bam O=$myname.gene_exon_tagged.bam ANNOTATIONS_FILE=$genomeDir/$genome.refFlat  TAG=GE
$dropseq/DetectBeadSynthesisErrors I=$myname.gene_exon_tagged.bam O=$myname.gene_exon_tagged.cleaned.bam OUTPUT_STATS=$myname.gene_exon_tagged.STATS SUMMARY=$myname.gene_exon_tagged.synthesis_stats.summary.txt NUM_BARCODES=$cellnumber2  PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
$dropseq/DigitalExpression I=$myname.gene_exon_tagged.cleaned.bam O=$myname.gene_exon_tagged.cleaned.dge.txt.gz SUMMARY=$myname.gene_exon_tagged.cleaned.summary.txt NUM_CORE_BARCODES=$cellnumber
$dropseq/BAMTagHistogram I=$myname.gene_exon_tagged.cleaned.bam O=$myname.gene_exon_tagged.cleaned.readscount.txt.gz TAG=XC
$mylib/Dosum.sh $myname
#less totalsummary | awk -F':' -v var2="$myname" -v var1="`date +%D`" -v var3="`tail -n +4 $myname.gene_exon_tagged.cleaned.summary.txt |awk '{if($3>1000) print $0}'|wc -l`" -v var4=`pwd` 'BEGIN{printf "%s\t%s\t",var1,var2} {printf "%i\t",$2 } END {printf "%i\t%s\n",var3,var4}' > $myname.dataQC.csv
Rscript --vanilla $mylib/QC.r $myname.gene_exon_tagged.cleaned.readscount.txt.gz  $cellnumber5 $myname.gene_exon_tagged.cleaned.summary.txt $myname.QC.pdf 1 2 help
#$mylib/ReadsSum.sh $myname ./
#for file in $myname.R1R2.sam $myname.R1R2.fastq $myname.unaligned_tagged_Cell.bam $myname.unaligned_tagged_CellMolecular.bam $myname.unaligned_tagged_filtered.bam $myname.unaligned_tagged_trimmed_smart.bam $myname.adapter_trimming_report.txt $myname.unaligned_mc_tagged_polyA_filtered.bam $myname.polyA_trimming_report.txt $myname.unaligned.header $myname.unaligned_mc_tagged_polyA_filtered.trimed.h.bam $myname.aligned.sorted.bam $myname.merged.bam $myname.gene_exon_tagged.STATS;do
#	rm $file
#done
#cp $myname.gene_exon_tagged.cleaned.dge.txt.gz ./DGE_files/
