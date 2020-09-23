#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyani

average_nucleotide_identity.py -o ../outputs/pyani_res -i ../input_genomes

