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

# This class provides some basic histogram reweighting functionality that can
# be used to obtain information from the macrostate distribution
class reweight(object):
    def __init__(self):
        self.set_num_smooth()
        self.set_phase_boundary()
        self.set_lnzs_guess()
        self._ignore_odd = False
        self._even_only = False
        self._saturation_found = False
        self._lnpi_vap_sat = None
        self._lnpi_liq_sat = None

    # find global minimum from this many macrostates
    def set_num_smooth(self, num_smooth=20):
        self._num_smooth = num_smooth

    # set the phase boundary
    def set_phase_boundary(self, phase_boundary=-1):
        if phase_boundary == -1:
            self._automatic_phase_boundary = True
            self._PHASE_BOUNDARY = -1
        else:
            self._automatic_phase_boundary = False
            self._PHASE_BOUNDARY = phase_boundary

    # get the last phase boundary
    def get_phase_boundary(self):
        return self._PHASE_BOUNDARY

    # set the guess for the activity(lnz) of saturation.
    # default: use the current activity
    def set_lnzs_guess(self, lnzsguess = -1):
        self._lnzsguess = lnzsguess

    # set if you would like to ignore odd macrostates
    def set_ignore_odd(self):
        self._ignore_odd = True

    # set if you have even macrostates only
    def set_even_only(self):
        self._even_only = True

    # obtain local minima with smoothing
    def _localMinSmooth(self, a):
        mins=[]
        for index in range(1, len(a)-1):
            #print index, len(a)
            lower = index-self._num_smooth
            if (lower<0): lower=0
            upper = index+self._num_smooth
            if (upper>= len(a)): upper = len(a)-1
            min = np.amin(a[lower:upper+1])
            if (math.fabs(min-a[index])<1e-7):
                #print index, min, lower, upper
                mins.append(index)
        return mins

    # return the number of particles which represents the phase boundary
    def _nboundary(self, lnpi):
        res=self._localMinSmooth(lnpi)
        #print(res)
        #assert(len(res)<2)
        if (len(res)==0):
            return None
        else:
            return res[0]

    # return the area
    def area(self, lnpi):
        sm=0
        for i in range(len(lnpi)):
            sm += math.exp(lnpi[i])
        return sm

    # return the normalized log Pi
    def norm(self, lnpi):
        lnpi = np.array(lnpi)
        if self._ignore_odd:
            import copy
            lnpi_working = [v for i, v in enumerate(lnpi) if i % 2 == 0]
        else:
            lnpi_working = lnpi
    #    print("shift", shift)
    #    lnpi -= shift
        shift = math.log(self.area(lnpi_working))
        lnpi -= shift
        return lnpi

    # separate the phase boundary and return the squared difference in the areas of these macrostate distributions
    def _sqprobdiff(self, criteria, lnzrw):
        #print("sqprw", lnzrw, lnz)
        criteria.lnPIrw(math.exp(lnzrw))
        lnpi = double_vector_to_list(criteria.lnPIrwdouble())
        if self._ignore_odd:
            lnpi = [v for i, v in enumerate(lnpi) if i % 2 == 0]
        if self._automatic_phase_boundary:
            pb=self._nboundary(lnpi)
        else:
            pb = self._PHASE_BOUNDARY
        # if no phase boundary, encourge multiple peaks by biasing toward similar
        # orders of magnitude of the first and last value of the lnpi
        if (pb == None):
            return (lnpi[0] - lnpi[-1])**2
        vap=lnpi[:pb]
        liq=lnpi[pb:]
        self._PHASE_BOUNDARY = pb
        self._lnpi_vap_sat = vap
        self._lnpi_liq_sat = liq
        #print("pb", pb)
        return (math.log(self.area(vap))-math.log(self.area(liq)))**2

    # find saturation by equating the probabilities of the two phases
    def lnzsat(self, criteria):
        from scipy.optimize import minimize
        lnzsguess = self._lnzsguess
        if lnzsguess == -1:
            lnzsguess = math.log(criteria.activ(0))
        res = minimize(lambda lnzrwi: self._sqprobdiff(criteria, lnzrw=lnzrwi[0]), lnzsguess, tol=1e-8)
        self._saturation_found = True
        return res["x"][-1]

    def average_macrostate(self, lnpi, first_macro=0):
        av = 0.
        fac = 1.
        if self._ignore_odd or self._even_only:
            fac = 2.
        for index, value in enumerate(lnpi):
            #print(index, value)
            av += fac*(index + first_macro)*math.exp(value)
        return av/self.area(lnpi)

    # return saturation properties
    def saturation_properties(self, space, criteria=None):
        if not self._saturation_found:
            assert(criteria != None) # provide criteria if saturation not already found
            lnzs = self.lnzsat(criteria)
            criteria.lnPIrw(math.exp(lnzs))
        lnpi = double_vector_to_list(criteria.lnPIrwdouble())
        return_dict = {"lnz_sat": math.log(criteria.activrw())}
        return_dict["phase_boundary"] = self.get_phase_boundary()
        return_dict["vapor"] = {"pressure": (-lnpi[0] + math.log(self.area(self._lnpi_vap_sat)))/space.volume()/criteria.beta(), "density": self.average_macrostate(self._lnpi_vap_sat)/space.volume()}
        return_dict["liquid"] = {"pressure": (-lnpi[0] + math.log(self.area(self._lnpi_liq_sat)))/space.volume()/criteria.beta(), "density": self.average_macrostate(self._lnpi_liq_sat, len(self._lnpi_vap_sat))/space.volume()}
        return return_dict
