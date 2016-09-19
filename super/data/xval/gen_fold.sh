#!/bin/sh

# e.g. ./gen_fold.sh promoters A [num_fold]

if [ $# -lt 2 ]
  then
    echo "Usage: ./gen_fold.sh <dataset> <split> [num_fold]"
    exit 1
fi

if [ $# -lt 3 ]
  then
    num_fold=10
  else
    num_fold=$3
fi

if [ ! -f ../raw/$1.in ]
  then
    echo "File not found: ../raw/$1.in"
    exit 1
fi

# create subdirectories, if necessary
mkdir -p $1
mkdir -p $1/$2

# make the necessary executables
cd ../../../src/base
make xval/gen_fold
make xval/cross_fold
make hercsearch
cd ../../super/data/xval/$1/$2

echo "generating $1.fold"
if [ ! -f $1$2.fold ]
  then
    ../../../../../bin/gen_fold $num_fold ../../../raw/$1.in > $1.fold
    rm -f $1${2}0.in
fi

echo "generating training files"
if [ ! -f $1${2}0.in ]
  then
    ../../../../../bin/cross_fold $num_fold ../../../raw/$1.in $1.fold
fi
