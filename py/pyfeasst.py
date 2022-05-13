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
    for _, vec1 in enumerate(vec):
        lst2 = list()
        for _, vec2 in enumerate(vec1):
          lst2.append(vec2)
        lst.append(lst2)
    return lst

def vector3d_to_list(vec):
    """ converts a swig stl vector to python list """
    lst = list()
    for _, vec1 in enumerate(vec):
        lst2 = list()
        for _, vec2 in enumerate(vec1):
            lst3 = list()
            for _, vec3 in enumerate(vec2):
                lst3.append(vec3)
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

# return equilibrium objective function
def equilibrium_objective(gce, beta_mu_rw, num_smooth):
    delta_conjugate = beta_mu_rw - gce.original_conjugate()
    gce.reweight(delta_conjugate)
    return gce.ln_prob().equilibrium_objective(num_smooth)

# find equilibrium
def find_equilibrium(gce, beta_mu_guess=-1, num_smooth=10):
    from scipy.optimize import minimize
    res = minimize(lambda beta_mu_rw: equilibrium_objective(gce, beta_mu_rw[0], num_smooth), beta_mu_guess, tol=1e-8)
    beta_mu_equilibrium = res["x"][-1]
    gce.reweight(beta_mu_equilibrium - gce.original_conjugate())
    return gce

# find the line in a csv file which begins with the characters 'state'
def line_beginning_with_state(csv_file):
    with open(csv_file, 'r') as f:
        lines = f.readlines()
        count = 0
        for line in lines:
            if line.split(',')[0] == 'state':
                return count
            count += 1
