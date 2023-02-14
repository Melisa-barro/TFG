#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=seq_toymodel
#SBATCH -t 00:30:00

module load gnu8/8.3.0

export HOME=/home/melisa/ACO_ADAPTADO

$HOME/aco $HOME/benchmarks/problema_99n.bs

