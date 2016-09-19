# runs a specified critic on a given image

here=$(dirname $0)

if [ $1 == "-c" ]
then
  f="-c"
  shift;
else
  f=""
fi

critic=$1
shift;
img=$1
shift;

metrics=$here/../bin/runMetrics

hercl=$here/../bin/hercl

if ! [ $critic ] || ! [ $img ]
then
  echo "Usage: run-critic.sh <critic.hrc> <img.png>"
  exit
fi

cmd="$hercl $* $critic +r <($metrics $f $img) <(echo )"
echo $cmd
eval $cmd
