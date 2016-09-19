#!/bin/sh

if [ $# -lt 1 ]
  then
    echo "Usage: ./run_many.sh <dataset> [target_cost]"
    exit 1
fi

if [ $# -lt 2 ]
  then
    cost=""
  else
    cost=$2
fi

for run in a b c d e f g h i j k l m n o
  do
    ../../run1.sh $1 $run $cost
done
