from pathlib import Path
import multiprocessing
import subprocess
import argparse
from pyfeasst import cd

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, help='FEASST directory')
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)

def run_file(filename):
    if 'checkpoint' not in filename.name:
        exclude = ['build', 'dev', 'feasst_test_env', 'library']
        if all([x not in str(filename.parent) for x in exclude]):
            with cd.cd(filename.parent):
                print("Running:", filename, "in", filename.parent)
                subprocess.call("rm tutorial_failures.txt tutorial_log.txt", shell=True, executable='/bin/bash')
                subprocess.call("jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=10000 --execute " + str(filename) + " > tutorial_log"+str(filename.name)+".txt 2>&1", shell=True, executable='/bin/bash')

def grep_file(filename):
    if 'checkpoint' not in filename.name:
        exclude = ['build', 'dev', 'feasst_test_env', 'library']
        if all([x not in str(filename.parent) for x in exclude]):
            with cd.cd(filename.parent):
                print("Grepping:", filename, "in", filename.parent)
                subprocess.call("grep \"Error\\|Assertion\" tutorial_log"+str(filename.name)+".txt >> tutorial_failures.txt", shell=True, executable='/bin/bash')
                subprocess.call("grep \"FAILED (fa\" " + str(filename) +" >> tutorial_failures.txt", shell=True, executable='/bin/bash')
                subprocess.call("grep \"Error\" " + str(filename) +" >> tutorial_failures.txt", shell=True, executable='/bin/bash')
                subprocess.call("grep \"ERROR\" " + str(filename) +" >> tutorial_failures.txt", shell=True, executable='/bin/bash')
                subprocess.call("grep \"feasst::CustomException\" " + str(filename) +" >> tutorial_failures.txt", shell=True, executable='/bin/bash')

#for filename in Path(feasst.install_dir()).rglob('*.ipynb'): run_file(filename)
pool = multiprocessing.Pool(4)
zip(pool.map(run_file, Path(ARGS.feasst_install).rglob('*.ipynb')))
pool = multiprocessing.Pool(1)
zip(pool.map(grep_file, Path(ARGS.feasst_install).rglob('*.ipynb')))
