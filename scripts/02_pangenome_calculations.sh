#!/bin/bash
set -e
# run this script with conda environment with gifrop install in it

# makes conda commands available in script 
# source ~/miniconda3/etc/profile.d/conda.sh

#activate conda environment with gifrop installed in it
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate gifrop

# make pangenomes directory
mkdir ./outputs/pan_genomes

# these two lines are to remove the header line from the metadata file
NUM_LINES=$(( $(cat ./outputs/01_metadata.tsv | wc -l) -1 ))
tail -n $NUM_LINES ./outputs/01_metadata.tsv > tempfile.tsv

while read a b c d e
  do

  # make pangenome directory and move appropriate genomes into each
  mkdir ./outputs/pan_genomes/"$e"
  cp ./outputs/annotations/"$b"/"$b".gff ./outputs/pan_genomes/"$e"/
  cp ./outputs/annotations/"$c"/"$c".gff ./outputs/pan_genomes/"$e"/
  cp ./outputs/annotations/"$d"/"$d".gff ./outputs/pan_genomes/"$e"/
  
  # enter the directory of the particular combination to calculate pangenome for
  cd ./outputs/pan_genomes/"$e"/
  
  # use gifrop package to:
  # 1) annotate with prokka
  # 2) calculate pangenome with roary:
  #       99% protein identity, output reference fasta containing all genes
  # 3) identify and extract 'genomic islands' with gifrop:
  #       min_genes = 1 so single genes can be islands
  #       also output separate fastas containing 1000bp of flanking DNA on either side of island
  
  # pan_pipe --roary_args '-p 20 -i 99 -e -n -z' --gifrop_args '--min_genes 1 --flank_dna 1000 --threads 20'
  roary -p 20 -i 99 -e -n -z *gff
  gifrop -m 1 --threads 20 --get_islands --flank_dna 1000
  # return to script directory
  cd -

done < tempfile.tsv


## ADD IN CALCULATION FOR ALL 6461s

mkdir ./outputs/pan_genomes/6461s

cp ./outputs/annotations/6461/6461.gff ./outputs/pan_genomes/6461s


while read a b c d e 
do
  if [ "$b" == '14229-5' -a "$c" == '6461' ]; then
    cp ./outputs/annotations/"$d"/"$d".gff ./outputs/pan_genomes/6461s/
  fi
done < tempfile.tsv


cd ./outputs/pan_genomes/6461s

roary -p 20 -i 99 -e -n -z *gff
gifrop -m 1 --threads 20 --get_islands --flank_dna 1000

abricate pan_genome_reference.fa > pan_genome_abricate.tsv
abricate --db vfdb pan_genome_reference.fa > pan_genome_vfdb.tsv

cd -


# ALL Genomes together in one big happy pangenome
mkdir ./outputs/pan_genomes/ALL

# find all gff annotation files and copy to ALL directory
find ./outputs/annotations/ -iname '*.gff' -exec cp {} ./outputs/pan_genomes/ALL \;

cd ./outputs/pan_genomes/ALL

roary -p 20 -i 99 -e -n -z *gff

abricate pan_genome_reference.fa > pan_genome_abricate.tsv
abricate --db vfdb pan_genome_reference.fa > pan_genome_vfdb.tsv

cd -

# 6461x13150s (and other tetO things)
mkdir ./outputs/pan_genomes/tetO

cp ./outputs/annotations/13150/13150.gff ./outputs/pan_genomes/tetO/
cp ./outputs/annotations/6461/6461.gff ./outputs/pan_genomes/tetO/

while read a b c d e 
do
  if [ "$b" == '6461' -a "$c" == '13150' ]; then
    cp ./outputs/annotations/"$d"/"$d".gff ./outputs/pan_genomes/tetO/
  fi
done < tempfile.tsv

# copy other 'inactive tetO' containing genomes in
find ./outputs/annotations/ -name 'GCA_*.gff' -exec cp {} ./outputs/pan_genomes/tetO/ \;

cd ./outputs/pan_genomes/tetO/

roary -p 20 -i 99 -e -n -z *gff

abricate pan_genome_reference.fa > pan_genome_abricate.tsv
abricate --db vfdb pan_genome_reference.fa > pan_genome_vfdb.tsv


cd -


rm tempfile.tsv

echo 'DONE!'
