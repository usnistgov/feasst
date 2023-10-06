from pathlib import Path
import multiprocessing
import subprocess
import feasst
from pyfeasst import cd

def run_file(filename):
    if 'checkpoint' not in filename.name:
        #if 'build' not in str(filename.parent):
        #exclude = ['build', 'dev', 'feasst_test_env']
        exclude = ['build', 'dev', 'feasst_test_env', 'library']
        if all([x not in str(filename.parent) for x in exclude]):
            with cd.cd(filename.parent):
                print("Running:", filename, "in", filename.parent)
                subprocess.call("rm tutorial_failures.txt tutorial_log.txt", shell=True, executable='/bin/bash')
                subprocess.call("jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=10000 --execute " + str(filename) + " > tutorial_log.txt 2>&1; grep \"Error\|Assertion\" tutorial_log.txt >> tutorial_failures.txt", shell=True, executable='/bin/bash')
                subprocess.call("grep \"FAILED (fa\" " + str(filename) +" >> tutorial_failures.txt", shell=True, executable='/bin/bash')
                subprocess.call("grep \"Error\" " + str(filename) +" >> tutorial_failures.txt", shell=True, executable='/bin/bash')
                subprocess.call("grep \"ERROR\" " + str(filename) +" >> tutorial_failures.txt", shell=True, executable='/bin/bash')
                subprocess.call("grep \"feasst::CustomException\" " + str(filename) +" >> tutorial_failures.txt", shell=True, executable='/bin/bash')

#for filename in Path(feasst.install_dir()).rglob('*.ipynb'): run_file(filename)
pool = multiprocessing.Pool(4)
zip(pool.map(run_file, Path(feasst.install_dir()).rglob('*.ipynb')))
