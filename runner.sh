#!/bin/bash
set -e

./00_setup.sh
Rscript 01_setup.R
./02_pangenome_calculations.sh
Rscript 03_HGT_ID.R
#knit report
Rscript 04_figures.R
Rscript 05_TetO_alignments.R
