######## Draft - Untested ########
#
##!/bin/sh
#
####################################################################
## pipeline for running read counting pipeline on 3' RNA-Seq data  #
##                                                                 #
## usage:  3prime_RNASeq_pipeline.sh  $basedir $num_of_processors  #
##                                                                 #
####################################################################
#
cd $1
CPU=$2
#
## bbduk trimming of raw fastq w/ fastqc before and after
#
## comment out for fixing STAR command
for file in `dir -d *.fastq.gz` ; do

    /home/jz322/miniconda3.7/bin/fastqc $file

    output=`echo "$file" | sed 's/R1.fastq.gz/trim.fastq.gz/'`

    /home/jz322/miniconda3.7/bin/bbduk.sh in=$file out=$output ftl=12 qtrim=rl trimq=20 minlen=35

    /home/jz322/miniconda3.7/bin/fastqc $output

done
#
## Prepare STAR genome:
## using gffread to convert gff3 to gtf format works!!!
/home/jz322/tools/STAR-2.7.5a/bin/Linux_x86_64_static/STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir ./  --genomeFastaFiles /home/jz322/Basil/3prime_RNASeq/Perrie_v1.0_chromosomes.fa --sjdbGTFfile /home/jz322/Basil/3prime_RNASeq/Perrie_braker_filter_gt.gtf --sjdbGTFtagExonParentTranscript Parent  --genomeSAindexNbases 13 --outSAMmapqUnique 60
#
## mapping using STAR:
## Note that running star with annotations is recommended, Star will extract splice junctions and improve alignment;
## In our 3' case, I don't think it hurts us anyway to do the alignment with or without - we tested this with pepper data and it didn't seem to make much difference
for file in `dir -d *trim.fastq.gz` ; do

    output2=`echo "$file" | sed 's/trim.fastq.gz/star.sort.bam/'`
    /home/jz322/tools/STAR-2.7.5a/bin/Linux_x86_64_static/STAR --runThreadN $CPU --outSAMtype BAM SortedByCoordinate --genomeDir ./  --sjdbGTFfile /home/jz322/Basil/3prime_RNASeq/Perrie_braker_filter_gt.gtf --sjdbGTFtagExonParentTranscript Parent  --readFilesIn $file --readFilesCommand zcat --outFileNamePrefix $output2 --outSAMmapqUnique 60

done
#
## Optional extension of existing features in annotation (if 3' ends aren't adequately annotated)
## This example extends each feature by 500 bp:
#
#bedtools slop -i Perrie_braker_filter_gt.gff -g sizes.genome -l 0 -r 500 -s > Perrie_braker_filter_gt_500.gff
#

# N.B. Htseq count is run here with -t gene, since gene features in annotations should include any 3' UTR features
#
for file in `dir -d *.bam` ; do
    outfile=`echo "$file" | sed 's/.sort.bamAligned.sortedByCoord.out.bam/.count/'`
    /home/jz322/miniconda3/bin/htseq-count -s yes -a 10 -t gene -i ID -f bam $file /home/jz322/Basil/3prime_RNASeq/Perrie_braker_filter_gt.gff > $outfile
done
