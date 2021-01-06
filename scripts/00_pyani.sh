#!/bin/bash

set -e 
echo 'extracting input genomes'
tar -xzvf input_genomes.tar.gz

parallel 'gunzip {}' ::: ./input_genomes/*gz
mkdir ./outputs/

echo 'calculating pairwise ANI values'
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyani

average_nucleotide_identity.py -o ./outputs/pyani_res -i ./input_genomes


conda deactivate

conda activate gifrop

mkdir ./outputs/abricate

# for x in ./input_genomes/*fna
#   do
#   sample=$(basename $x .fna)
#   abricate $x > ./outputs/abricate/"$sample".abricate
#   abricate --db vfdb $x > ./outputs/abricate/"$sample".vfdb
# done
echo 'annotating input genomes with prokka'
# AMR ID
parallel 'abricate {} > ./outputs/abricate/{/.}.abricate' ::: ./input_genomes/*fna
# Virulence ID
parallel 'abricate --db vfdb {} > ./outputs/abricate/{/.}.vfdb' ::: ./input_genomes/*fna


#annotate with prokka
parallel 'prokka -outdir ./outputs/annotations/{/.} -prefix {/.} --proteins ./C_jejuni_NCTC_11168.gbff {}' ::: ./input_genomes/*fna
