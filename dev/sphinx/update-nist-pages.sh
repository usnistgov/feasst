#!/bin/bash

# For preparing a new release:
# Check /feasst/dev/analyze_public_interface.py
# CMakeLists.txt:
# - version
# - compiler flags (-Wall -pedantic -g)
# - remove depend.py
# add feasst.h
# tag commit
# run nightly test

# update nist-pages
mkdir build
cd build
cmake -DUSE_SPHINX=ON ..
make html > tt 2>&1
grep -v "_arguments.rst: WARNING: document" tt | grep -v "_arguments.rst:4: WARNING: Duplicate"  | grep -v "^Declaration is" | grep -v "WARNING: Duplicate C++ declaration, also defined"
version=$(git describe)
#branch=`git branch | grep \* | cut -d ' ' -f2`
mv html html2
git checkout nist-pages
cp -r html2/* ../
git status
echo "Press [Enter] to commit with message: $version?"
read -rs
git commit -a -m "$version"
#git checkout $branch
