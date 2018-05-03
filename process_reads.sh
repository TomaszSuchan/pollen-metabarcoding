#!/usr/bin/env bash
# the raw reads should be in directory called "1-rawreads"

mkdir 2-merged/
mkdir 3-cleaned/
mkdir 4-filtered/
mkdir 5-dereplicated/
mkdir reads-classified_sintax
THREADS=4 # change according to the number of threads available on your machine

#merge reads
for fvd in 1-rawreads/*_R1_001.fastq.gz
do
rev=${fvd/R1/R2}
out=$(basename $fvd | sed 's/_R1_001.fastq.gz//g')
pearRM -j $THREADS -f $fvd -r $rev -o 2-merged/$out  2>&1 | tee 2-merged/$out.log
done

#remove primer sequences
for f in 2-merged/*.assembled.fastq
do
out=$(basename $f| sed 's/.assembled.fastq//g')
read_fastq -e base_33 -i $f | \
remove_primers -f ATGCGATACTTGGTGTGAAT -r GCATATCAATAAGCGGAGGA -m 0 -i 0 -d 0 | \
write_fastq -x -o 3-cleaned/$out.cleaned.fastq
done

#filter reads
for f in 3-cleaned/*.cleaned.fastq
do
out=$(basename $f| sed 's/.cleaned.fastq//g')
vsearch --threads $THREADS \
        --fastq_filter  $f \
        --fastq_maxee 0.5 \
        --fastq_minlen 250 \
        --fastq_maxns 0 \
        --fastaout 4-filtered/$out.filtered.fasta \
        --fasta_width 0 2>&1 | tee 4-filtered/$out.filtered.log
done

#dereplicate, discard reads which appear once, rename reads with sample id:
for f in 4-filtered/*.filtered.fasta
do
out=$(basename $f| sed 's/.filtered.fasta//g')
label=${out/-/}
vsearch --threads $THREADS \
        --derep_fulllength $f \
        --minuniquesize 2 \
        --strand plus \
        --output 5-dereplicated/$out.derep.fasta \
        --sizeout \
        --uc 5-dereplicated/$out.derep.uc \
        --relabel $label. \
        --fasta_width 0 2>&1 | tee 5-dereplicated/$out.derep.log
done

#SINTAX classification
for f in 5-dereplicated/*.derep.fasta
do
out=$(basename $f .derep.fasta)
usearch10 -sintax $f \
          -db viridiplantae_all_2014.sintax.fa \
          -tabbedout reads-classified_sintax/$out.sintax.txt \
          -strand both \
          -sintax_cutoff 0.95
done