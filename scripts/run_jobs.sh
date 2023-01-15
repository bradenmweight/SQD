#!/bin/bash

for n in {5..199}; do 
    cp NAMD.in TRAJ/traj-${n}/; 
    cp submit.MMMM_NAMD TRAJ/traj-${n}/; 
    cd TRAJ/traj-${n}/; 
        sbatch submit.MMMM_NAMD ; 
        cd ../../; done