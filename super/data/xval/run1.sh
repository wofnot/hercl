#!/bin/sh

# e.g. ./run1.sh promoters 15

if [ $# -lt 2 ]
  then
    echo "Usage: ./run1.sh <dataset> <num_runs>"
    exit 1
fi

DIR=$PWD
cd ../../../../../src/base
make hercsearch
cd $DIR

mkdir -p c0$2
for fold in 0 1 2 3 4 5 6 7 8 9
do
  if [ -f $1$fold.in ]
  then
    echo "../../../../../bin/hercsearch -s $1$fold.in -n $cost > c0$2/$1$fold.out"
    ../../../../../bin/hercsearch -s $1$fold.in -n $cost > c0$2/$1$fold.out
  fi
done

echo "Done."
