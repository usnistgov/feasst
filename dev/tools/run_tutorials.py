from pathlib import Path
import multiprocessing
import feasst
import pyfeasst

def run_file(filename):
    if 'checkpoint' not in filename.name:
        if 'build' not in str(filename.parent):
            with pyfeasst.cd(filename.parent):
                print("Running:", filename, "in", filename.parent)
                pyfeasst.bash_command("rm tutorial_failures.txt")
                pyfeasst.bash_command("jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=10000 --execute " + str(filename))
                pyfeasst.bash_command("grep \"FAILED (fa\" " + str(filename) +" >> tutorial_failures.txt")
                pyfeasst.bash_command("grep \"Error\" " + str(filename) +" >> tutorial_failures.txt")

#for filename in Path(feasst.install_dir()).rglob('*.ipynb'): run_file(filename)
pool = multiprocessing.Pool(4)
zip(pool.map(run_file, Path(feasst.install_dir()).rglob('*.ipynb')))
