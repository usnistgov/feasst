from pathlib import Path
import multiprocessing
import subprocess
from pyfeasst import cd

for filename in Path('../').rglob('launch*.py'):
    if 'checkpoint' not in filename.name:
        if 'build' and 'dev' and 'feasst_test_env' not in str(filename.parent):
            with cd.cd(filename.parent):
                print("Running:", filename.name, "in", filename.parent)
                subprocess.call("python " + str(filename.name) + " 2> launch.log", shell=True, executable='/bin/bash')
#                subprocess.call("jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=10000 --execute " + str(filename) + " > tutorial_log.txt 2>&1; grep \"Error\|Assertion\" tutorial_log.txt >> tutorial_failures.txt")
#                subprocess.call("grep \"FAILED (fa\" " + str(filename) +" >> tutorial_failures.txt")
#                subprocess.call("grep \"Error\" " + str(filename) +" >> tutorial_failures.txt")
#                subprocess.call("grep \"ERROR\" " + str(filename) +" >> tutorial_failures.txt")
#                subprocess.call("grep \"feasst::CustomException\" " + str(filename) +" >> tutorial_failures.txt")

# put all ids in nums list
nums=list()
for filename in Path('../').rglob('launch_ids.txt'):
    print(filename, filename.name, filename.parent)
    with open(filename) as f:
        nums.append(f.read().splitlines())
#print(nums)

# search squeue for ids. If no ids are present, then jobs are finished. Proceed
finished = False
while not finished:
    subprocess.call("squeue -u $LOGNAME > squeue.txt", shell=True, executable='/bin/bash')
    with open('squeue.txt') as f: squeue = f.read().splitlines()
    #print(squeue)
    finished = True
    for nms in nums:
        for num in nms:
            for line in squeue:
                if num in line:
                    finished = False
                    print('num', num, 'line', line)
            # if num in squeue.txt
    import time
    time.sleep(60)

# scrap job output for errors
for filename in Path('../').rglob('hostname_*.out'):
    with cd.cd(filename.parent):
        subprocess.call("grep \"Err\" " + str(filename.name) + " >> launch_failures.txt", shell=True, executable='/bin/bash')
        subprocess.call("grep \"Throw\" " + str(filename.name) + " >> launch_failures.txt", shell=True, executable='/bin/bash')
        subprocess.call("grep \"No such file or directory\" " + str(filename.name) + " >> launch_failures.txt", shell=True, executable='/bin/bash')
    subprocess.call("grep \"Err\" " + str(filename) + " >> launch_failures.txt", shell=True, executable='/bin/bash')
    subprocess.call("grep \"Throw\" " + str(filename) + " >> launch_failures.txt", shell=True, executable='/bin/bash')
    subprocess.call("grep \"No such file or directory\" " + str(filename) + " >> launch_failures.txt", shell=True, executable='/bin/bash')

for filename in Path('../').rglob('launch.log'):
    with cd.cd(filename.parent):
        subprocess.call("grep \"Err\" " + str(filename.name) + " >> launch_failures.txt", shell=True, executable='/bin/bash')
        subprocess.call("grep \"FAILED\" " + str(filename.name) + " >> launch_failures.txt", shell=True, executable='/bin/bash')
    subprocess.call("grep \"Err\" " + str(filename) + " >> launch_failures.txt", shell=True, executable='/bin/bash')

