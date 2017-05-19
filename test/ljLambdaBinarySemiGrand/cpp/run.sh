prog=semigrand
FEASST_INSTALL_DIR_=$HOME/jeetain/feasstv0.1.3/
echo $FEASST_INSTALL_DIR_
$FEASST_INSTALL_DIR_/tools/compile.sh $prog
mkdir -p tmp #directory for checkpoint files
./$prog
