prog=muvttmmclj
~/feasst/tools/compile.sh $prog
cp ~/feasst/src/main $prog
export OMP_NUM_THREADS=4
mkdir -p tmp #directory for checkpoint files
./$prog
