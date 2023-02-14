#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=seq_toymodel
#SBATCH --gres=gpu
#SBATCH -t 00:30:00

module load gnu8/8.3.0
module load cuda/10.2.89

export HOME=/home/melisa/ACO_GPU_2

echo CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES

$HOME/aco $HOME/benchmarks/problema_4155n.bs

