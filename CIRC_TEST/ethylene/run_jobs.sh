#!/bin/bash

for n in {0..100}; do 
    cp NAMD.in TRAJ/traj-${n}/; 
    cp submit.SQD TRAJ/traj-${n}/; 
    cd TRAJ/traj-${n}/; 
        rm -rf MD* output.slurm #G16 EL_STRUCTURE*
        sbatch submit.SQD ; 
        cd ../../; done
