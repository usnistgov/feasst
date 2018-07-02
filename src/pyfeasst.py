"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

# PyFeasst - collection of python utility functions
import os, sys
#from libxdrfile import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc,
#from libxdrfile import DIM, exdrOK
import feasst
import numpy as np

def read_log_file(file_name, xcol=0, ycol=1):
    """ Read log file from already-run analysis
    Skip commented lines beginning with "#"
    file_name:   name of log file
    xcol:       column for x data
    ycol:       column for y data"""
    xlst = list()
    ylst = list()
    log = open(file_name)
    lines = log.readlines()
    for i in range(len(lines)):
        columns = lines[i].split(' ')
        if (columns[0] != "#"):
            xlst.append(columns[xcol].rstrip())
            ylst.append(columns[ycol].rstrip())
    return xlst, ylst

def string_array_to_np(xarr, yarr):
    """ convert string array, often obtained by read_log_file,
    into a numpy array of floats """
    tmp = [float(val) for val in xarr]
    xdata = np.array(tmp)
    tmp = [float(val) for val in yarr]
    ydata = np.array(tmp)
    return xdata, ydata

def log2np(file_name, xcol=0, ycol=1):
    """ Interface to convert output of read_log_file
    into numpy arrays of floats
    file_name:   name of log file
    xcol:       column for x data
    ycol:       column for y data"""
    xlst, ylst = read_log_file(file_name, xcol, ycol)
    return string_array_to_np(xlst, ylst)

def plot(name, xlst, ylst, xlabel="", ylabel="", xlabelfs=18, ylabelfs=28):
    """ Plot x,y data into name.eps,jpg """
    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt
    #fig = plt.figure(figsize=size)
    plt.plot(xlst, ylst, 'k-')
    plt.plot(xlst, ylst, 'ko')
    plt.xlabel(xlabel, fontsize=xlabelfs)
    plt.ylabel(ylabel, fontsize=ylabelfs)
    plt.tight_layout()
    plt.savefig(name+".eps")
    plt.savefig(name+".jpg")

def list_to_double_vector(lst):
    """ converts a python list to a swig stl vector """
    vec = feasst.DoubleVector(len(lst))
    for index in range(len(lst)):
        vec[index] = float(lst[index])
    return vec

def list_to_int_vector(lst):
    """ converts a python list to a swig stl vector """
    vec = feasst.IntVector(len(lst))
    for index in range(len(lst)):
        vec[index] = int(lst[index])
    return vec

def double_vector_to_list(vec):
    """ converts a swig stl vector to python list """
    lst = list()
    for index in range(len(vec)):
        lst.append(vec[index])
    return lst

class cd:
    """Context manager for changing the current working directory
       Thanks to Brian Hunt on Stack Overflow
       http://stackoverflow.com/questions/431684/how-do-i-cd-in-python/24176022
       #24176022 """
    def __init__(self, newPath):
        self._newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self._savedPath = os.getcwd()
        os.chdir(self._newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self._savedPath)

import math
def string2sumSq(text):
    rtrn = 0
    for val in text.split():
        rtrn += float(val)**2
    return rtrn
assert(math.fabs(1 - string2sumSq("0.46187167792334405 -0.35043322571889868 -0.22710692686778397 0.78249188571715755")) < feasst.DTOL)
