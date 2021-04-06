#!/bin/bash
set -e

./scripts/00_pyani.sh
Rscript ./scripts/01_setup.R
./scripts/02_pangenome_calculations.sh
Rscript ./scripts/03_HGT_ID.R
#knit report
Rscript ./scripts/04_figures.R
Rscript ./scripts/05_TetO_alignments.R
Rscript ./scripts/06_Resistance_annotations.R
Rscript ./scripts/07_rpsL_SNPs_compare.R

