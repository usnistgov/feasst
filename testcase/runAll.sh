smryfile=`pwd`/summary.log
rm $smryfile
pushd lj/srsw/nvt-mc
  ./run.sh >> $smryfile
popd


