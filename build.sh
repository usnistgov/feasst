mkdir -p build
cd build
cmake -DUSE_SWIG=ON ..
make -j 4
