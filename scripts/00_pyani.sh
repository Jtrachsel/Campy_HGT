#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyani

average_nucleotide_identity.py -o ../outputs/pyani_res -i ../input_genomes


conda deactivate

conda activate gifrop

mkdir ./outputs/abricate

for x in ./input_genomes/*fna
  do
  sample=$(basename $x .fna)
  abricate $x > ./outputs/abricate/"$sample".abricate
done
