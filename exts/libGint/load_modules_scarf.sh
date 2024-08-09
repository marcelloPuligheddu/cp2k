srun --partition=gpu-devel --gres=gpu:1 --mem=0 --ntasks-per-node=1 -C "scarf21" --time=8:00:00  --pty /bin/bash

module use /work4/scd/scarf562/eb-common/modules/all
module load amd-modules OpenMPI/4.1.4-GCC-11.3.0 NCCL/2.12.12-GCCcore-11.3.0-CUDA-11.7.1
module load UCX-CUDA/1.12.1-GCCcore-11.3.0-CUDA-11.7.1 ScaLAPACK/2.2.0-gompi-2022a-fb NVHPC/23.1-CUDA-12.0.0

