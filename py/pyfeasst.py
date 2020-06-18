# pyfeasst - collection of python utilities

import os
#import feasst

class cd:
    """Context manager for changing the current working directory
       Thanks to Brian Hunt on Stack Overflow
       https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python/13197763#13197763 """
    def __init__(self, newPath):
        self._newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self._savedPath = os.getcwd()
        os.chdir(self._newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self._savedPath)

def vector2d_to_list(vec):
    """ converts a swig stl vector to python list """
    lst = list()
    for index1 in range(len(vec)):
        lst2 = list()
        for index2 in range(len(vec[index1])):
          lst2.append(vec[index1][index2])
        lst.append(lst2)
    return lst

def vector3d_to_list(vec):
    """ converts a swig stl vector to python list """
    lst = list()
    for index1 in range(len(vec)):
        lst2 = list()
        for index2 in range(len(vec[index1])):
            lst3 = list()
            for index3 in range(len(vec[index1][index2])):
                lst3.append(vec[index1][index2][index3])
            lst2.append(lst3)
        lst.append(lst2)
    return lst

import subprocess
def bash_command(cmd):
    """Execute a bash command using subprocess"""
    subprocess.call(cmd, shell=True, executable='/bin/bash')
    #subprocess.Popen(['/bin/bash', '-c', cmd])

def read_checkpoint(filename):
    """Return contents of checkpoint file as a string"""
    with open (filename, "r") as myfile:
        checkpoint=myfile.readlines()
    assert(len(checkpoint) == 1)  # checkpoint files should have only one line
    return checkpoint[0]

#def forcefield_dir(filename=''):
#    return os.path.dirname(os.path.realpath(feasst.__file__)) + '/forcefield/' + filename
