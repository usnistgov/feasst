prog=lj
$FEASST_INSTALL_DIR_/tools/compile.sh $prog
cp $FEASST_INSTALL_DIR_/src/main $prog
mkdir -p tmp #directory for checkpoint files
./$prog
