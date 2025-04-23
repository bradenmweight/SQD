#!/bin/bash
#for NMOL in 1 2 5 10 25 50; do
for NMOL in 5 10 25 50; do
    sbatch submit.spectra $NMOL 
done