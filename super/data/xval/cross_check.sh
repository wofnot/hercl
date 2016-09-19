#!/bin/sh

if [ $# -lt 3 ]
  then
    echo "Usage: ./cross_check <dataset> <num_fold> <num_runs>"
    exit 1
fi

DIR=$PWD
cd ../../../../../src/base
make xval/cross_check
cd $DIR

../../../../../bin/cross_check ../../../raw/$1.in $2 $1.fold $3 ${1}*.hrc

