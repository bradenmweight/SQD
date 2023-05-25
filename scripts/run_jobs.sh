#!/bin/bash

for n in {0..50}; do 
    cp NAMD.in TRAJ/traj-${n}/; 
    cp submit.SQD TRAJ/traj-${n}/; 
    cd TRAJ/traj-${n}/; 
        rm -rf MD* output.slurm
        sbatch submit.SQD ; 
        cd ../../; done