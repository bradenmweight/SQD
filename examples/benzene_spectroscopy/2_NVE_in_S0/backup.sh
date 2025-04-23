#!/bin/bash

for n in {0..200}; do
    echo "Backing up TRAJ/traj-${n}/"
    cd TRAJ/traj-${n}/
    git add . --all --force
    cd ../../
    sleep 0.05
done