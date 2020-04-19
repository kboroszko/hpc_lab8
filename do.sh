#!/bin/bash
make clean
git pull
make
sbatch generate.batch
squeue -u kb385496
