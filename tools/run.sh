prog=nvt
mkdir -p tmp
$FEASST_INSTALL_DIR_/tools/compile.sh $prog
cp $FEASST_INSTALL_DIR_/src/main $prog
export OMP_NUM_THREADS=1
./$prog > $prog.log 

