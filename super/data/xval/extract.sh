#!/bin/sh

for run in a b c d e f g h i j k l m n o
do
  for fold in 0 1 2 3 4 5 6 7 8 9
  do
    if [ -f c0$run/*$fold.out ]
      then
        echo "=A$fold$run"
        tail -3 c0$run/*$fold.out | head -2
    fi
  done
done
