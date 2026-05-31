#!/bin/sh
###############################################################
# pipeline for running RNAseq DE on multiple paired end files #
#                                                             #
# usage:                                                      #
#                                                             #
#      RNAseq_pe_pipeline.sh $indir $CPU $bin                 #
#                                                             #
###############################################################
#get commandline options
cd /home/jz322/Basil/trichomes_RNAseq/RNAseq_data    #input file
CPU=20   #no of CPUS
bin=/home/jz322/chloroseq_se_pipeline   #directory containing executables
#ref=$4   #path and name of reference fasta
#gtf=$5   #path and name of reference gtf

#filter out rRNA
#for file in `dir -d *_R1.fastq` ; do
#    file2=`echo "$file" | sed 's/_R1/_R2/'`
#    merge=`echo "$file" |sed 's/.fastq/.merged.fastq/'`
#    $bin/sortmerna-2.1b/scripts/merge-paired-reads.sh $file $file2 $merge
#done

#ls *merged.fastq |parallel -j 1 $bin/sortmerna-2.1b/sortmerna --ref $bin/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,$bin/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95-db:$bin/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,$bin/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98-db --reads {} --aligned {.}_rRNA --other {.}_other --num_alignments 1 --fastx --paired_in -a $CPU --log

#for file in `dir -d *other.fastq` ; do
#    read1=`echo "$file" | sed 's/other/R1_other/'`
#    read2=`echo "$file" | sed 's/other/R2_other/'`
#    $bin/sortmerna-2.1b/scripts/unmerge-paired-reads.sh $file $read1 $read2
#done

#clean reads
#for file in `dir -d *_R1_other.fastq`; do
#   file2=`echo "$file" | sed 's/_R1_/_R2_/'`
#   out1=`echo "$file" |sed 's/.fastq/_cln.fastq/'`
#   unpaired1=`echo "$file" |sed 's/.fastq/_cln_unpaired.fastq/'`
#   out2=`echo "$file2" |sed 's/.fastq/_cln.fastq/'`
#   unpaired2=`echo "$file2" |sed 's/.fastq/_cln_unpaired.fastq/'`

#   java -jar $bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE $file $file2 $out1 $unpaired1 $out2 $unpaired2 ILLUMINACLIP:$bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:50

#done

#get read count 
#nohup wc *cln.fastq | awk '{print $1/4 "\t" $4}' > filtered_read_count.txt 2>&1 & 

#map reads to genome
for file in `dir -d *_1.fastq` ; do

    #change "*R1*.fastq_" to "*R1*clean.fq"                                                                                          
    file2=`echo "$file" |sed 's/_1/_2/'`
    samfile=`echo "$file" | sed 's/_1.fastq/.sam/'`
    
    hisat2 --max-intronlen 120000 --dta -p $CPU -x /home/jz322/Basil/trichomes_RNAseq/Perrie -1 $file -2 $file2 -S $samfile

done

##convert sam files and get stats
ls *.sam |parallel --gnu -j $CPU samtools view -Sb -o {.}.bam {}
ls *.bam |parallel --gnu -j $CPU samtools sort -o {.}.sort.bam {}
ls *.sort.bam |parallel --gnu -j $CPU samtools flagstat {} ">" {.}.flagstat
cat *flagstat |grep "mapped (" |sed 's/.*(\(.*\)%.*/\1/g' |awk '{sum+=$1} END { print "Average = ",sum/NR}' > average_mapping.txt

##calc mapping efficiency
ls *sort.bam |parallel -j $CPU --pipe samtools view -F 260 {} |wc ">" mapped_reads.txt

##calc mapping percent
for file in `dir -d *.flagstat` ; do

    num_mapped=`echo "$file" | sed 's/.flagstat/.mapped/'`
    grep "mapped (" $file | sed 's/.*(\(.*\)%.*/\1/g' > $num_mapped

done

for file in `ls -1 *mapped`; do
    echo "$file" > ./tmpfile
    cat "$file" >> ./tmpfile
    mv ./tmpfile "$file"
    cat $file |tr "\n" "\t" |sed -e '$a\' > $file.csv
done
cat *.csv > mapping_efficiency.csv

##get fpkm and counts
#mkdir ballgown

for file in `dir -d *.sort.bam` ; do
    
    outfile=`echo "$file" | sed 's/.bam/.gtf/'`  
    outdir=`echo "$file" |sed 's/.bam//'`
    
    stringtie  -e -B -p $CPU -G /home/jz322/Basil/3prime_RNASeq/Perrie_braker_filter_gt.gtf -o ballgown/$outdir/$outfile $file

done

python3 $bin/prepDE.py -i ballgown -g gene_count_matrix.csv -t transcript_count_matrix.csv 

##import to R with DESeqDataSetFromMatrix (DEseq) and DGEList (EdgeR)
