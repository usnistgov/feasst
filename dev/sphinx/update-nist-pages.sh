#!/bin/bash

# For preparing a new release:
# bump version number
# Check /feasst/dev/tools/analyze_public_interface.py
# CMakeLists.txt:
# - check version, check -DEV=ON, /feasst/dev/tools/depend.py
# update html and nist-pages
# run nightly test
# tag commit

# update nist-pages
mkdir build
cd build
cmake -DUSE_SPHINX=ON ..
make html > tt 2>&1
grep -v "_arguments.rst: WARNING: document" tt | grep -v "_arguments.rst:4: WARNING: Duplicate"  | grep -v "^Declaration is" | grep -v "WARNING: Duplicate C++ declaration, also defined" | grep -v "/home/hwh/feasst/README.rst: WARNING: document isn't included in any toctree" | grep -v "warning: The following parameters of feasst::" | grep -v "warning: The following parameter of feasst::"
version=$(git describe)
#branch=`git branch | grep \* | cut -d ' ' -f2`
cp -r html/* ../html/
git commit
cd ..
cp -r html html2
git checkout nist-pages
cp -r html2/* .
git status
echo "Press [Enter] to commit with message: $version?"
read -rs
git commit -a -m "$version"
#git checkout $branch
