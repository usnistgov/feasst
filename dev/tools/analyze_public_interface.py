import json
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--output_interface", "-p", help="file name to output interface json file", type=str, default='public_interface.json')
parser.add_argument("--compare_previous_interface", "-c", help="json file of p", type=str, default="")
args = parser.parse_args()
print('args', args)

def analyze_header_file(filename):
    """ Obtain public interface arguments from header, assuming one public class per header. """
    first_public_arguments = "  /** @name Arguments"
    last_public_arguments = "  //@}"
    with open(filename, 'r') as fl:
        lines = fl.readlines()
    num = 0
    targs = list()
    classname = ""
    for line in lines:
        if num == 1:
            #print(line[:6])
            if line[:6] == "    - ":
                ln = line[6:]
                if ":" in ln:
                    ln = ln[:ln.index(":")]
                    #print('ln', ln)
                    targs.append(ln)
                elif " arguments" in ln:
                    if " arguments." in ln:
                        ln = ln[:ln.index(" arguments.")]
                        #print('ln', ln)
                        targs.append(ln)
                    else:
                        ln = ln[:ln.index(" arguments")]
                        #print('ln', ln)
                        targs.append(ln)
                else:
                    #print('ln', ln)
                    assert False, "unrecognized"
                    targs.append(ln)
        if first_public_arguments in line:
            num += 1
            #print(filename)
        if last_public_arguments in line:
            break
        #print(line[:6])
        if line[:6] == "class ":
            ln = line[6:]
            if " " in ln:
                ln = ln[:ln.index(" ")]
                #print('ln', ln)
                classname = ln
        if num > 1:
            print("Multiple public interface arguments in the same header file:", filename)
            exit()
    #print('classname', classname, 'targs', targs)
    return classname, targs

# obtain the public interface as a dictionary of headers and arguments
public = {}
print('\n')
print('***********************************')
print('classes containing public arguments')
print('***********************************')
print('\n')
for filename in Path('../plugin/').rglob('*.h'):
    if 'build' and 'dev' and 'feasst_test_env' and 'library' not in str(filename.parent):
        #print(filename)
        classname, targs = analyze_header_file(filename)
        if len(targs) > 0:
            print(classname, filename)
            public[classname] = targs
#print(public)
with open(args.output_interface, 'w') as fl:
    fl.write(json.dumps(public))

# search for nonunique classname
for index1, class1 in enumerate(public):
    for index2, class2 in enumerate(public):
        if index1 != index2 and class1 == class2:
            print(class1, class2)
            assert False, "nonunique class name shouldn't be possible."

# search for nonunique arguments (OK)
for class1 in public:
    for class2 in public:
        if class1 != class2:
            #print(class1, class2)
            for arg in public[class1]:
                if arg in public[class2]:
                    if arg not in public:
                        assert True
                        print('shared argument', arg, 'in', class1, class2)

import dictdiffer
if args.compare_previous_interface != "":
    with open(args.compare_previous_interface, 'r') as fl:
        previous = json.load(fl)
    #print('previous', previous)
    print('\n')
    print('***********')
    print('differences')
    print('***********')
    print('\n')
    for diff in list(dictdiffer.diff(previous, public)):
        print(diff)
