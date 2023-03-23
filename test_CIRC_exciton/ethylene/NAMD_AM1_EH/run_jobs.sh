#!/bin/bash

# Set something like {0..200} to {0..500} for Braden-scale submissions
for n in {0..1}; do 
    cp NAMD.in TRAJ/traj-${n}/ 
    cp submit.SQD TRAJ/traj-${n}/ 
    cd TRAJ/traj-${n}/ 
        rm -rf G16 MD* output.slurm
        #sbatch submit.SQD 
        cd ../../
done