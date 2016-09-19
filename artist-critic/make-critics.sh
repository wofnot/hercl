# generates a bunch of critics for the given phase

folder=$1
shift;

phase=$1
num_critics=$2
tolerance=$3

if ! [ $phase ] || ! [ $folder ]
then
  echo "Usage: make-critics.sh <folder> <phase #> [<num_critics>]"
  exit
fi

if ! [ $num_critics ]
then
  num_critics=1
fi

if ! [ $tolerance ]
then
  tolerance=0.1
fi



critics=$folder/$phase/critics

mkdir -p $critics

hercsearch=$folder/../../bin/hercsearch

if ! [ -f $hercsearch ]
then
  echo "couldn't find hercsearch executable at $hercsearch"
  exit 1;
fi

flags="-n -t $tolerance"

data=$folder/$phase/train.in

for (( i=1; i<=$num_critics; i++ ))
do
  if ! [ -f $critics/$i.hrc ]
  then
    # run hercsearch and extract solution from output
    echo generating critic $i
    $hercsearch $flags -d $i -s $data | tee /dev/tty | sed -e 's/^TI/`TI/' -e 's/^\*\*\*/`\*\*\*/' | awk '/`TI/{p=1;print;next} p&&/^Done/{p=0};p' > $critics/$i.hrc
  else 
    echo $critics/$i.hrc already exists, skipping
  fi
done

