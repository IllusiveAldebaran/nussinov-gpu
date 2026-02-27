#!/bin/sh

repoDir=$(dirname "$(realpath "$0")")
echo $repoDir
cd $repoDir

make clean
make run