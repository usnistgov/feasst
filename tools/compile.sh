driver=$1.cc
source=$FEASST_INSTALL_DIR_/src
cp $driver $source/main.cc
cd $source
rm main
make main
