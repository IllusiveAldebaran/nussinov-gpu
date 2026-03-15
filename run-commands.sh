#!/bin/sh

repoDir=$(dirname "$(realpath "$0")")
echo $repoDir
cd $repoDir

make clean
#make run CPU=1 SEQ_N=25

# CPU run
make CPU=1 SEQ_N=1500
./nussinov synthetic/synthetic_sequences.txt > cpu_seq1500.txt

#make clean
#
## GPU run (number of threads parameter)
#make BLOCK=512
#nsys profile --stats=true --force-overwrite=true ./nussinov $(head -c 1500 inputs/ec16s.seq) > gpu_thread512.txt
