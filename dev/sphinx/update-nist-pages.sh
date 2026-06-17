#!/bin/bash

# For preparing a new release:
# bump version number
#  - located both in CMakeLists.txt and pyproject.toml
# test branch using nightly build and test / launch
# Check /feasst/dev/tools/analyze_public_interface.py
# CMakeLists.txt:
# - check version, check -DEV=ON, /feasst/dev/tools/depend.py
# Use cmake -DUSE_PIP=ON -DEV=ON .. to:
#       /feasst/dev/tools/inheritance.py
#       /feasst/feasst python menu.py --particles ../build/public_particles.txt --descripts ../build/public_description.json --factory ../build/public_factories.json ->> updates /feasst/feasst/data/menu.json file
# update html and nist-pages
# tag commit

# update nist-pages
mkdir build
cd build
cmake -DUSE_SPHINX=ON ..
make html > tt 2>&1
grep -v "_arguments.rst: WARNING: document" tt | grep -v "_arguments.rst:4: WARNING: Duplicate"  | grep -v "^Declaration is" | grep -v "WARNING: Duplicate C++ declaration, also defined" | grep -v "/home/$LOGNAME/feasst/README.rst: WARNING: document isn't included in any toctree" | grep -v "warning: The following parameters of feasst::" | grep -v "warning: The following parameter of feasst::" | grep -v "_skbuild"
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

## update pypi https://packaging.python.org/en/latest/tutorials/packaging-projects/
# cd ~/feasst
# rm -r _skbuild
# rm -r src/feasst.get-info
##python3 -m pip install --upgrade build
# Ensure there is no MANIFEST.in that is messing up files
#CMAKE_BUILD_PARALLEL_LEVEL=8 python3 -m build --sdist # source-only distribution (or make wheels using cibuildwheel)
##CMAKE_BUILD_PARALLEL_LEVEL=8 CMAKE_ARGS="-DUSE_PYBIND11=ON" python3 -m build
#python3 -m twine check dist/*
#python3 -m pip install --upgrade twine rich
#python3 -m twine upload --repository testpypi dist/*
#
##open a fresh install with wsl
#cmd
#wsl --unregister Ubuntu-26.04 # remove old install
#wsl --install Ubuntu-26.04
#cd
#sudo apt update
#sudo apt upgrade
## test build instructions in documentation, but install from pypi
#CMAKE_BUILD_PARALLEL_LEVEL=8 pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple feasst==0.25.19.4
#feasst
#feasst-menu

