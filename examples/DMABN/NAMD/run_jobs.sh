#!/bin/bash

for n in {201..400}; do 
    cp NAMD.in TRAJ/traj-${n}/; 
    cp submit.MMMM_NAMD TRAJ/traj-${n}/; 
    cd TRAJ/traj-${n}/; 
        rm -rf G16 MD* output.slurm
        sbatch submit.MMMM_NAMD ; 
        cd ../../; done