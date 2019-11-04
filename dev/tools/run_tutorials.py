from pathlib import Path
import feasst
import pyfeasst

pyfeasst.bash_command("rm tutorial_failures.txt")
for filename in Path(feasst.install_dir()).rglob('*.ipynb'):
    if 'checkpoint' not in filename.name:
      with pyfeasst.cd(filename.parent):
          print("Running:", filename, "in", filename.parent)
          pyfeasst.bash_command("jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=10000 --execute " + str(filename))
          pyfeasst.bash_command("grep \"FAILED (fa\" " + str(filename) +" >> tutorial_failures.txt")
          pyfeasst.bash_command("grep \"Error\" " + str(filename) +" >> tutorial_failures.txt")
