"""
Given a FEASST text input file, search the documentation to output possibly
relevant references for the methods used.
"""

import argparse
import depend

def parse():
    """ Parse arguments """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_file', type=str, help='FEASST text input file', required=True)
    parser.add_argument('--source_dir', type=str, help='/path/to/feasst', required=True)
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    print(args)
    return args

def used_classes(filename):
    """
    Read all lines in input file, remove all empty and comment lines
    and obtain the first space-separated value of each line as classes.
    """
    with open(str(filename), 'r') as file1:
        lines = file1.readlines()
    cls = list()
    for _, line in reversed(list(enumerate(lines))):
        if line[0] == '#' or line[0] == '\n':
            lines.remove(line)
        else:
            cls.append(line.split(' ')[0].strip('\n'))
    return cls

def find_class_refs(args, depends):
    class_refs = list()
    for include in depends:
        #print(include[0])
        header_file_name = args.source_dir + '/plugin/' + include[0]
        cls = depend.read_class(header_file_name)
        with open(header_file_name) as header_file:
            lines = header_file.readlines()
        refs = list()
        for line in lines:
            if ':footcite:' in line:
                start = line.find(':footcite:')
                cutline = line[(start+13):]
                #print('cutline:', cutline)
                end = cutline.find('`')
                #print('end', end)
                cutline = cutline[:end]
                refs.append(cutline)
        if len(refs) > 0:
            class_refs.append([cls, refs])
    return class_refs

def read_references(args):
    """ Read bibtex and return entries by name """
    bibpath = args.source_dir + '/dev/sphinx/refs.bib'
    with open(bibpath, 'r') as bibfile:
        lines = bibfile.readlines()
    names = list()
    entries = list()
    refs = {}
    current_name = ''
    entry = ''
    for line in lines:
        new_entry = False
        for beginning in ['article', 'book']:
            cut = len(beginning)+2
            if line[:cut] == '@'+beginning+'{':
                new_entry = True
                new_cut = cut
        if new_entry:
            current_name = line[new_cut:-2]
            names.append(current_name)
            if len(entry) > 1:
                refs[previous_name] = entry
                entry = ''
            previous_name = current_name
        entry += line
    refs[previous_name] = entry
    return refs

def user_print(used_classes, references_by_class, references):
    """ Print the references for classes used in input file. """
    print('# Recommended citation(s) for FEASST')
    print(references['hatch_monte_2024'])
    #print('used_classes', used_classes)
    unique_classes = list(set(used_classes))
    #print('unique_classes', unique_classes)
    for used in unique_classes:
    #for used in used_references_by_class:
        #print('used', used)
        for index, cls in enumerate(references_by_class):
            #print('checking cls', cls, ':', cls[0])
            if used in cls[0]:
                print('# Recommended citation(s) for', used)
                for ref in references_by_class[index][1]:
                    if ref in references:
                        print(references[ref])
                    else:
                        print('ref:', ref, 'is missing from /feasst/dev/sphinx/refs.bib')
        #    print(references[ref])

if __name__ == '__main__':
    ARGS = parse()
    CLS = used_classes(ARGS.input_file)
    INCLUDE_PLUGIN = depend.read_plugins(ARGS.source_dir+'/build/plugins.txt')
    DEPS = depend.dependency(ARGS, INCLUDE_PLUGIN)
    CLS_REF = find_class_refs(ARGS, DEPS)
    print(CLS_REF)
    REFS = read_references(ARGS)
    user_print(CLS, CLS_REF, REFS)
    #print('refs', REFS)
    #print('000')
    #print(REFS['metropolis_equation_1953'])
    #print('001')
    #print(REFS['errington_direct_2003'])
    #print(REFS['frenkel_understanding_2002'])
