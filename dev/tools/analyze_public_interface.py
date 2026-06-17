import json
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--output_interface", "-p", help="file name to output interface json file", type=str, default='public_interface.json')
parser.add_argument("--output_descript", "-d", help="file name to output description json file", type=str, default='public_description.json')
parser.add_argument("--output_particles", "-o", help="file name to output particle file names", type=str, default='public_particles.txt')
parser.add_argument("--compare_previous_interface", "-c", help="json file of p", type=str, default="")
args = parser.parse_args()
print('args', args)

def analyze_header_file(filename):
    """ Obtain public interface arguments from header, assuming one public class per header. """
    #print('reading', filename)
    first_public_arguments = "  /** @name Arguments"
    last_public_arguments = "  //@}"
    with open(filename, 'r') as fl:
        lines = fl.readlines()
    num = 0
    is_descript_found = False
    targs = list()
    targdes = list()
    classname = ""
    descript = ""
    is_arg_des_found = False
    for line in lines:
        if num == 1:
            #print(line[:6])
            if line[:6] == "    - ":
                ln = line[6:]
                if ":" in ln:
                    idx = ln.index(":")
                    l = ln[:idx]
                    #print('l', l)
                    targs.append(l)
                    targdes.append(ln[(idx+2):])
                    is_arg_des_found = True
                elif " arguments" in ln:
                    if " arguments." in ln:
                        ln = ln[:ln.index(" arguments.")]
                        #print('ln', ln)
                        targs.append(ln)
                        targdes.append("")
                        is_arg_des_found = True
                    else:
                        ln = ln[:ln.index(" arguments")]
                        #print('ln', ln)
                        targs.append(ln)
                        targdes.append("")
                        is_arg_des_found = True
                else:
                    #print('ln', ln)
                    assert False, "unrecognized"
            elif is_arg_des_found:
                if line[:5] == "   */":
                    is_arg_des_found = False
                else:
                    #print('line', line, len(targdes))
                    #exit()
                    targdes[-1] += line
        if first_public_arguments in line:
            num += 1
            #print(filename)
        if last_public_arguments in line:
            break
        if line[:3] == " */":
            is_descript_found = False
        if line[:3] == "/**":
            is_descript_found = True
        elif is_descript_found:
            descript += line
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
    def href_to_http(descript, parent):
        # determine directories back
        parent = str(parent)[3:]
        nback=0
        if descript.count('a href') == 1:
            nback = descript.count('../')
        elif descript.count('a href') > 1:
            assert False, descript
        start = '<a href="'
        end='.html">'
        start_idx = descript.find(start) + len(start)
        end_idx = descript.find(end, start_idx)
        newrs = ""
        newr = ''
        result=""
        if start_idx != -1 and end_idx != -1:
            result = descript[start_idx:end_idx]
            paths = str(parent).split('/')[:-nback]
            for pth in paths:
                newrs += pth+'/'
            paths = result.split('/')[nback:]
            for pth in paths:
                newr += pth+'/'
            newr = newr[:-1]
            descript = descript.replace(result, newrs+newr)
        #descript += ':parent:'+str(parent)+':nback:'+str(nback)+':newrs:'+str(newrs)+':result:'+str(result)+":newr:"+str(newr)
        descript = descript.replace('<a href="', '(https://pages.nist.gov/feasst/').replace('.html">', '.html)').replace('</a>', '')
        return descript
    descript = descript.replace('\n\n', '__TEMP1234__').replace('\n', ' ').replace('__TEMP1234__', '\n\n').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('\n ', '\n')
    descript = href_to_http(descript, filename.parent)
    if descript != "":
        if descript[0] == ' ':
            descript = descript[1:]
    #print('classname', classname, 'targs', targs)
    for it,tad in enumerate(targdes):
        tad = href_to_http(tad, filename.parent)
        #tad = tad.replace('<a href=', '(https://pages.nist.gov/feasst/').replace('</a>', ')')
        targdes[it] = tad.replace('\n', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ')
    return classname, targs, descript, targdes

# obtain the public interface as a dictionary of headers and arguments
public = {}
doc = {}
print('\n')
print('***********************************')
print('classes containing public arguments')
print('***********************************')
print('\n')
for filename in Path('../plugin/').rglob('*.h'):
    #if 'build' and 'dev' and 'feasst_test_env' and 'library' not in str(filename.parent):
    excluded = ['html', 'build', 'dev', 'feasst_test_env', 'library', '_sources', '_skbuild']
    if not any(word in str(filename.parent) for word in excluded):
        #print(filename)
        classname, targs, descript, targdes = analyze_header_file(filename)
        if len(targs) > 0:
            print(classname, filename)
            public[classname] = targs
        # separate required arguments from optinal arguments, based if the description has a key word
        assert len(targs) == len(targdes)
        required_args = []
        optional_args = []
        for idx in range(len(targs)):
            if targdes[idx] != '' and 'optional' not in targdes[idx] and 'Optional' not in targdes[idx] and 'default' not in targdes[idx]:
                required_args.append([targs[idx], targdes[idx]])
            else:
                optional_args.append([targs[idx], targdes[idx]])
        doc[classname] = {'descript': descript, 'required': required_args, 'optional': optional_args}
#print(public)
with open(args.output_interface, 'w') as fl:
    fl.write(json.dumps(public))
with open(args.output_descript, 'w') as fl:
    fl.write(json.dumps(doc, indent=2))

# search for nonunique classname
for index1, class1 in enumerate(public):
    for index2, class2 in enumerate(public):
        if index1 != index2 and class1 == class2:
            print(class1, class2)
            assert False, "nonunique class name shouldn't be possible."

## search for nonunique arguments (OK)
#for class1 in public:
#    for class2 in public:
#        if class1 != class2:
#            #print(class1, class2)
#            for arg in public[class1]:
#                if arg in public[class2]:
#                    if arg not in public:
#                        assert True
#                        print('shared argument', arg, 'in', class1, class2)

if args.compare_previous_interface != "":
    import dictdiffer
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

# output list of particle files
particle_files = dict()
for filename in Path('../').rglob('particle/*.txt'):
    excluded = ['html', 'build', 'dev', 'feasst_test_env', 'library', '_sources', '_skbuild', 'dist']
    if not any(word in str(filename.parent) for word in excluded):
    #if 'html' and 'build' and 'dev' and 'feasst_test_env' and 'library' and '_sources' and '_skbuild' not in str(filename.parent):
        #print(filename, filename.parent)
        #print(str(filename)[3:])
        with open(filename, 'r') as file1:
            text = ''
            descript = file1.read().splitlines()
            for line in descript:
                text += line.replace('\n\n', '__TEMP1234__').replace('\n', ' ').replace('__TEMP1234__', '\n\n').replace('\n ', '\n') + '\n'
        particle_files[filename.name] = [text, str(filename)[3:]]
        #particle_files.append(str(filename)[3:])

#particle_files = sorted(particle_files)
#particle_files.sort(key=str.lower)
particle_files = {k: particle_files[k] for k in sorted(particle_files)}
frequent = ['atom_new.txt', 'lj_new.txt', 'spce_new.txt', 'tip4p.txt', 'hard_sphere.txt', 'dimer_mie_CO2.txt', 'dimer_mie_N2.txt', 'co2.txt', 'ethane.txt', 'atom2d.txt', 'janus.txt', 'mie.txt', 'mie_C1H4_methane.txt', 'mie_C2H6_ethane.txt' ]
for freq in reversed(frequent):
    if freq in particle_files:
        des = particle_files[freq]
        particle_files.pop(freq)
        particle_files = {freq: des, **particle_files}

with open(args.output_particles, 'w') as fl:
    fl.write(json.dumps(particle_files))
    #for line in particle_files:
    #    fl.write(line+'\n')

#print('sorted', particle_files)
