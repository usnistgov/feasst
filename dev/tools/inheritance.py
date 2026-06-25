import json
from pathlib import Path

def string_between_strings(line, begin, end):
    """ Return the string that is between beginning and end strings """
    return line[line.find(begin)+len(begin):line.find(end)]

# search all headers for derived and base class tree
depend = {}
factories = []
first = 'class '
middle = ' : public '
last1 = ' {\n'
last2 = ',\n'
first_fac = '  std::shared_ptr<'
last_fac = '> factory(const std::string name'
for filename in Path('../plugin/').rglob('*.h'):
    if 'build' and 'dev' and 'feasst_test_env' and 'library' not in str(filename.parent):
        #print(filename)
        with open(filename, 'r') as fl:
            lines = fl.readlines()
        for iline,line in enumerate(lines):
            # find inheritance
            if line[:len(first)] == first and middle in line:
                if line[-len(last1):] == last1:
                    #print(line)
                    derived=string_between_strings(line, first, middle)
                    base=string_between_strings(line, middle, last1)
                    assert ' ' not in derived and ' ' not in base, base+'->'+derived
                    depend[derived]=[base]
                if line[-len(last2):] == last2:
                    if lines[iline+1][-len(last2):] == last2:
                        assert False, "triple inheritance not implemented."
                    elif lines[iline+1][-len(last1):] == last1:
                        #print(line)
                        derived=string_between_strings(line, first, middle)
                        base1=string_between_strings(line, middle, last2)
                        base2=string_between_strings(lines[iline+1], 'public ', last1)
                        assert ' ' not in derived and ' ' not in base1 and ' ' not in base2, base1+','+base2+'->'+derived
                        depend[derived]=[base1,base2]
                    else:
                        assert False, 'unrecognized:'+lines[iline+1]
            # find factories
            if line[:len(first_fac)] == first_fac and last_fac in line:
                factories.append(string_between_strings(line, first_fac, last_fac))
                #print(factories[-1])
                #print(line)#, filename)

def add_to_tree(tree, branch, leaf):
    if branch in tree:
        tree[branch].append(leaf)
    else:
        tree[branch] = [leaf]

base_tree = {}
for fac in factories:
    for dep in depend:
        if fac in depend[dep]:
            add_to_tree(base_tree, fac, dep)
            for dep2 in depend:
                if dep in depend[dep2]:
                    add_to_tree(base_tree, fac, dep2)
                    for dep3 in depend:
                        if dep2 in depend[dep3]:
                            add_to_tree(base_tree, fac, dep3)
                            for dep4 in depend:
                                if dep3 in depend[dep4]:
                                    add_to_tree(base_tree, fac, dep4)
                                    for dep5 in depend:
                                        if dep4 in depend[dep5]:
                                            add_to_tree(base_tree, fac, dep5)
                                            for dep6 in depend:
                                                if dep5 in depend[dep6]:
                                                    assert False, 'implement further inheritance depth'

# alphabetize
for fac in base_tree:
    base_tree[fac].sort()

# Metropolis first in Criteria factory, AlwaysReject last
crit = base_tree['Criteria']
if 'Metropolis' in crit:
    crit.insert(0, crit.pop(crit.index('Metropolis')))
if 'AlwaysReject' in crit:
    crit.append(crit.pop(crit.index('AlwaysReject')))

#print(base_tree['Model'], len(base_tree['Model']))
#print(depend)
with open('public_factories.json', 'w') as file1:
    file1.write(json.dumps(base_tree))
