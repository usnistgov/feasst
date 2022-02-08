from pathlib import Path
import multiprocessing
#import feasst
import pyfeasst

for filename in Path('../').rglob('launch*.py'):
    if 'checkpoint' not in filename.name:
        if 'build' and 'dev' and 'feasst_test_env' not in str(filename.parent):
            with pyfeasst.cd(filename.parent):
                print("Running:", filename.name, "in", filename.parent)
                pyfeasst.bash_command("python " + str(filename.name) + " 2> launch.log")
#                pyfeasst.bash_command("jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=10000 --execute " + str(filename) + " > tutorial_log.txt 2>&1; grep \"Error\|Assertion\" tutorial_log.txt >> tutorial_failures.txt")
#                pyfeasst.bash_command("grep \"FAILED (fa\" " + str(filename) +" >> tutorial_failures.txt")
#                pyfeasst.bash_command("grep \"Error\" " + str(filename) +" >> tutorial_failures.txt")
#                pyfeasst.bash_command("grep \"ERROR\" " + str(filename) +" >> tutorial_failures.txt")
#                pyfeasst.bash_command("grep \"feasst::CustomException\" " + str(filename) +" >> tutorial_failures.txt")

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
    pyfeasst.bash_command("squeue -u $LOGNAME > squeue.txt")
    with open('squeue.txt') as f: squeue = f.read().splitlines()
    #print(squeue)
    finished = True
    for num in nums[0]:
        for line in squeue:
            if num in line:
                finished = False
                print('num', num, 'line', line)
        # if num in squeue.txt
    import time
    time.sleep(60)

# scrap job output for errors
for filename in Path('../').rglob('hostname_*.out'):
    with pyfeasst.cd(filename.parent):
        pyfeasst.bash_command("grep \"Err\" " + str(filename.name) + " >> launch_failures.txt")
    pyfeasst.bash_command("grep \"Err\" " + str(filename) + " >> launch_failures.txt")

for filename in Path('../').rglob('launch.log'):
    with pyfeasst.cd(filename.parent):
        pyfeasst.bash_command("grep \"Err\" " + str(filename.name) + " >> launch_failures.txt")
        pyfeasst.bash_command("grep \"FAILED\" " + str(filename.name) + " >> launch_failures.txt")
    pyfeasst.bash_command("grep \"Err\" " + str(filename) + " >> launch_failures.txt")

