#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=seq_toymodel
#SBATCH -t 00:30:00

module load gnu8/8.3.0

export HOME=/home/patricia/ACO_GPU/toymodel_aco_seq/aco/

$HOME/aco $HOME/benchmarks/problema_4155n.bs

