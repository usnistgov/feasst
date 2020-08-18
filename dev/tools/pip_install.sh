#!/bin/bash

set +x
# activate a new env
python3 -m venv ~/feasst_env
source ~/feasst_env/bin/activate  # can leave by deactivate

# clean env?
# rm -r ~/feasst_env/lib/python3.7/site-packages/feasst-*.egg/

rm -r build
mkdir -p build
cd build

ln -s ../forcefield .
cmake -DUSE_SWIG=ON ..
make _feasst -j$CPU_COUNT

# here , you may want to update the version in setup.py
#cd ..
#ln -s build feasst

pip install wheel twine
# readme issues? pip install readme_renderer
#python3 build/setup.py sdist bdist_wheel
python3 setup.py sdist bdist_wheel
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

# install into env
pip install -i https://test.pypi.org/simple/ feasst==0.8.6

# cleanup

rm -r feasst dist feasst.egg-info
