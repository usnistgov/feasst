# Custom script to launch jobs
#!/bin/bash
ppn=12
nodes=1
((procs=$nodes*$ppn))
if [ $ppn == 12 ]; then
  ppnvis=24
else
  ppnvis=$ppn
fi
output="$(basename $0 .sh)"

for temp in 1.5; do
  case $temp in
    1.5) x=370; lnz=-1.568214; p=0;
  esac
  prog=muvt

  dir=t${temp}x${x}z${lnz}p$p
  mkdir -p $dir
  cd $dir
  if [ `hostname` == "raritan.nist.gov" ]; then
    dir=`pwd | sed 's/home/wrk/'`
    #queue="-q medium"
    queue="-q long"
    #queue="-q small"
    nodeexclus=
    resource="nodes=${nodes}:ppn=$ppnvis:rp12"
    resource="nodes=${nodes}:ppn=$ppnvis:ivy30"
  else
    dir=`pwd`
    nodeexclus="-n"
    resource="nodes=${nodes}:ppn=$ppnvis"
  fi
  mkdir -p $dir
  output2=${output}
  $FEASST_INSTALL_DIR_/tools/compile.sh ../$prog
  cp ../$prog.cc $dir
  cp $FEASST_INSTALL_DIR_/src/main.cc $dir
  cp $FEASST_INSTALL_DIR_/src/main $dir/$prog

# write output script to file
cat << _EOF_ > ${output2}.cmd
#!/bin/bash
#PBS -l $resource
#PBS $nodeexclus
#PBS $queue
#PBS -V
mkdir -p $dir/
cd $dir
mkdir -p tmp
echo "Running on host \$(hostname)"
echo "Time is \$(date)"
echo "Directory is \$PWD"

##move to tmp directory and run
#TMPDIR=/tmp/$LOGNAME/\$PBS_JOBID
#echo "TMPDIR is \$TMPDIR"
#mkdir -p \$TMPDIR
#cp -r $dir/* \$TMPDIR
#cd \$TMPDIR
#~/tools/rsyncslepper.sh . $dir/ &

export OMP_NUM_THREADS=$ppn
./$prog -x $x -z $lnz > ${prog}.log

##copy files to wrk
#cp -r * $dir/
#rm -r *

echo "Job is done"
echo "Time is \$(date)"
_EOF_

  jobid=`qsub ${output2}.cmd`
  cd ..
done
