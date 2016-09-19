# generates the training data for specified phase, accumulating negative and positive items from all previous phases
folder=$1
shift;

phase=$1
shift;

use_colour=$1
shift;

extra_folders=$*
if ! [ $phase ] || ! [ $folder ]
then
  echo "Usage: make-train.sh <folder> <phase #>"
  exit
fi

if [ $use_colour -eq 0 ]
then
  flags=""
else
  flags="-c"
fi

mkdir -p $folder/$phase

metrics=$folder/../../bin/runMetrics

data=$folder/$phase/train.in

pos=$folder/$phase/pos.tmp
neg=$folder/$phase/neg.tmp

# accumulate all positive items from previous phases into pos.in
> $pos
for (( x=1; x<=$phase; x++ ))
do
  $metrics $flags $folder/$x/pos/* >> $pos
done

# accumulate all negative items
> $neg
for (( x=1; x<=$phase; x++ ))
do
  $metrics $flags $folder/$x/neg/* >>  $neg
done

# also grab items from any extra folders

for extra in $extra_folders
do
  $metrics $flags $extra/neg/* >> $neg
  $metrics $flags $extra/pos/* >> $pos
done

n_pos=`wc -l < $pos`
n_neg=`wc -l < $neg`

echo "writing to $data";

> $data

# pos items
for (( x=1; x<=$n_pos; x++ ))
do
  head -$x $pos | tail -1 >> $data
  echo ,1 >> $data
done
# neg items
for (( x=1; x<=$n_neg; x++ ))
do
  head -$x $neg | tail -1 >> $data
  echo ,0 >> $data
done

if [ $n_pos -gt $n_neg ]
then
  src=$neg
  dif=$((n_pos-n_neg))
  val=0
  n_src=$n_neg
else
  src=$pos
  dif=$((n_neg-n_pos))
  val=1
  n_src=$n_pos
fi

# 1 -> 2
# 98 -> 99
# 99 -> 100
# 100 -> 1

# extra items from whichever one needs more
for (( x=1; dif>0; x=(x%$n_src + 1) ))
do
  ((dif--))
  
  head -$x $src | tail -1 >> $data
  echo ,$val >> $data
done

# OLD METHOd:
## put random selection of 300 pos and 300 neg items into train.in

#> $data
#for (( x=1; x<=$items; x++ ))
#do
#  shuf -n 1 $pos >> $data
#  echo ,1 >> $data
#  shuf -n 1 $neg >> $data
#  echo ,0 >> $data
#done

