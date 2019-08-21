'''
  This utility automatically runs all tutorials
'''

import os
import argparse
parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')
required.add_argument("--source_dir", "-s", help="/path/to/feasst", type=str, required=True)
args = parser.parse_args()
import pyfeasst


for dir_, _, files in os.walk(args.source_dir + '/plugin/'):
    # run cpp files
    if 'main.cpp' in files:
        print(dir_, files)
        with pyfeasst.cd(dir_):
            os.mkdir('build')
            with pyfeasst.cd('build'):
                pyfeasst.bash_command('cmake .. && make && ./main')

    # run py files
    pys = [k for k in files if 'py' in k]
    for pyfile in pys:
        print('pyfile', pyfile, dir_)
        with pyfeasst.cd(dir_):
            pyfeasst.bash_command('../../../py/run.sh ' + pyfile)
