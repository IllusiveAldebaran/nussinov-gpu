#!/bin/sh

repoDir=$(dirname "$(realpath "$0")")
echo $repoDir
cd $repoDir

make clean
make run CPU=1 SEQ_N=25