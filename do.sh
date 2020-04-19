#!/bin/bash
make clean
git pull
make
sbatch graph.batch
squeue -u kb385496
