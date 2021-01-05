#!/bin/bash
set -e

source ~/README_FUNCTIONS

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/945/185/GCA_004945185.1_PDT000199611.2/GCA_004945185.1_PDT000199611.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/005/266/575/GCA_005266575.1_PDT000199856.2/GCA_005266575.1_PDT000199856.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/947/465/GCA_004947465.1_PDT000199954.2/GCA_004947465.1_PDT000199954.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/947/585/GCA_004947585.1_PDT000199985.2/GCA_004947585.1_PDT000199985.2_genomic.fna.gz

gunzip *gz

for x in *fna; do mv $x "${x%%.*}".fna; done

parallel 'rename_contigs {} {.}' ::: *fna

mv *fna ../input_genomes/

