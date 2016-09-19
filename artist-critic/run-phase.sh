here=$(dirname $0)

folder=$1
shift;
phase=$1

if ! [ $phase ] || ! [ $folder ]
then
  echo "Usage: run-phase.sh <folder> <phase #>"
  exit
fi

# read num_artists from settings.txt
settings=$folder/settings.txt

artist_tolerance=`grep artist_tolerance $settings | cut -d':' -f2`
push_pop=`grep push_pop $settings | cut -d':' -f2`
no_loop=`grep no_loop $settings | cut -d':' -f2`
max_step=`grep max_step $settings | cut -d':' -f2`
max_movement=`grep max_movement $settings | cut -d':' -f2`
use_colour=`grep use_colour $settings | cut -d':' -f2`

# write art-school settings to art-school.txt
as_settings=$folder/art-school.txt
> $as_settings
echo "L critics.list" >> $as_settings
echo "out.png" >> $as_settings
echo $artist_tolerance >> $as_settings
echo $push_pop >> $as_settings
echo $no_loop >> $as_settings
echo $max_step >> $as_settings
echo $max_movement >> $as_settings
echo $use_colour >> $as_settings

critic_tolerance=`grep critic_tolerance $settings | cut -d':' -f2`
num_artists=`grep num_artists $settings | cut -d':' -f2`
num_critics=`grep num_critics $settings | cut -d':' -f2`
use_past_critics=`grep use_past_critics $settings | cut -d':' -f2`
n_cells=`grep cells $settings | cut -d':' -f2`
train_data=`grep train_data $settings | cut -d':' -f2`
max_eval=`grep max_eval $settings | cut -d':' -f2`

echo "running phase $folder:$phase."
for cmd in "$here/make-train.sh $folder $phase $use_colour $train_data" "$here/make-critics.sh $folder $phase $num_critics $critic_tolerance" "$here/make-artists.sh $folder $phase $num_artists $use_past_critics $n_cells $max_eval"
do
  echo $cmd
  $cmd
done
