#!/bin/sh

repoDir=$(dirname "$(realpath "$0")")
echo $repoDir
cd $repoDir

make clean

# CPU & GPU run
make CPU=1 SEQ_N=800 BLOCK=512 GRID=8
nsys profile --stats=true ./nussinov synthetic/synthetic_sequences.txt > results.txt

# Extract CPU lines
sed '/Running again on GPU/q' results.txt | grep "^('" > cpu_only.txt
# Extract GPU lines
sed -n '/Running again on GPU/,$p' results.txt | grep "^('" > gpu_only.txt

# Compare 
if diff cpu_only.txt gpu_only.txt > /dev/null; then
    echo "---------------------------------------"
    echo "CPU and GPU outputs match."
    echo "---------------------------------------"
else
    echo "---------------------------------------"
    echo "CPU and GPU outputs mismatch."
    echo "---------------------------------------"
    
    diff -y --suppress-common-lines cpu_only.txt gpu_only.txt | head -n 10
fi


