# generates a bunch of artists for the given phase
# also saves their output to the neg directory for the next phase

folder=$1
shift;

phase=$1
num_artists=$2
use_past_critics=$3
n_cells=$4
max_eval=$5

if ! [ $phase ] || ! [ $folder ] 
then
  echo "Usage: make-artists.sh <folder> <phase #> [num_artists] [use_past_critics] [n_cells] [max_eval]"
  exit
fi

if ! [ $num_artists ]
then
  num_artists=1
fi
if ! [ $use_past_critics ]
then
  use_past_critics=0
fi
if ! [ $n_cells ]
then
  n_cells=1
fi

if ! [ $max_eval ]
then
  max_eval=99000000000
fi

artists=$folder/$phase/artists

mkdir -p $artists

artschool=../../bin/art-school # relative to $folder
debug_artschool=../../bin/debug_art-school # relative to $folder
if ! [ -f $folder/$artschool ]
then
  echo "couldn't find art-school executable at $artschool"
  exit 1;
fi

flags="-n -k 1 -z 0 -c $n_cells -g $max_eval"

# set phase in source code

## critic=$phase/critics/1.hrc # relative to $folder

settings=$folder/art-school.txt

if ! [ -f $settings ]
then
  echo "couldn't find settings file at $settings"
  exit 1;
fi

critic_list=$phase/critics.list # relative to $folder

if [ $use_past_critics -eq 0 ]
then 
  e=$phase
else
  e=1
fi

# add past and new critics to critic_list
cd $folder
> $critic_list
w=1 # weighting
for (( p=$phase; p>=$e; p-- ))
do
  n_critics=`ls $p/critics | wc -l`
  for critic in `ls $p/critics/*`
  do
    t=`echo "scale=10; $w/$n_critics" | bc`
    echo $t $critic >> $critic_list
  done 
  w=`echo "scale=10; $w/2" | bc`
done

cd -

# update art-school.txt to use critic list
sed -i x '1s:.*:L '$critic_list':' $settings

lib=$phase/artists.lib

cd $folder
> $lib
# generate library of all past artists, relative to $folder
for (( p=1; p<$phase; p++ ))
do
  ls $p/artists/* >> $lib
done
cd -

output_folder=$((phase+1))/neg # relative to $folder
mkdir -p $folder/$output_folder
for (( i=1; i<=$num_artists; i++ ))
do
  if ! [ -f $artists/$i.hrc ]
  then
    
    fn=$output_folder/$i.png
    sed -i x '2s:.*:'$fn':' $settings
     
    cd $folder
    
    # run artschool and extract solution from output
    echo generating artist $i
    # log=`$artschool $flags -d $i -L $lib > /dev/tty`
    
    artist_file=./$phase/artists/$i.hrc
    
    # TI if finished, MI if max_eval exceeded
    echo $artschool $flags -d $i -L $lib
    $artschool $flags -d $i -L $lib | tee /dev/tty | sed -e 's/^\([TM]I\)/`\1/' -e 's/^\*\*\*/`\*\*\*/' | awk '/`[TM]I/{p=1;print;next} p&&/^Done/{p=0};p' > $artist_file
    
    if grep MI $artist_file
    then
      # max_eval exceeded so we need to render the image'
      echo rendering
      $debug_artschool $artist_file -b
    fi
    
    # jump back
    cd -
  else 
    echo $artists/$i.hrc already exists, skipping
  fi
done

