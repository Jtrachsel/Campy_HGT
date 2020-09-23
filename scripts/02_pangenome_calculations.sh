#!/bin/bash

# run this script with conda environment with gifrop install in it

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gifrop

# make pangenomes directory
mkdir ../outputs/pan_genomes

# these two lines are to remove the header line from the metadata file
NUM_LINES=$(( $(cat ../outputs/01_metadata.tsv | wc -l) -1 ))
tail -n $NUM_LINES ../outputs/01_metadata.tsv > tempfile.tsv

while read a b c d 
  do

  # make pangenome directory and move appropriate genomes into each
  mkdir ../outputs/pan_genomes/"$d"
  cp ../input_genomes/"$a".fna ../outputs/pan_genomes/"$d"/
  cp ../input_genomes/"$b".fna ../outputs/pan_genomes/"$d"/
  cp ../input_genomes/"$c".fna ../outputs/pan_genomes/"$d"/
  
  # enter the directory of the particular combination to calculate pangenome for
  cd ../outputs/pan_genomes/"$d"/
  
  # use gifrop package to:
  # 1) annotate with prokka
  # 2) calculate pangenome with roary:
  #       99% protein identity, output reference fasta containing all genes
  # 3) identify and extract 'genomic islands' with gifrop:
  #       min_genes = 1 so single genes can be islands
  #       also output separate fastas containing 1000bp of flanking DNA on either side of island
  
  pan_pipe --roary_args '-p 20 -i 99 -e -n -z' --gifrop_args '--min_genes 1 --flank_dna 1000 --threads 20'
  
  # return to script directory
  cd -

done < tempfile.tsv

rm tempfile.tsv

