import numpy as np
from libxdrfile import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc, DIM, exdrOK

xtc = '1L2Y.xtc'

# get number of atoms
natoms = read_xtc_natoms(xtc)
print natoms

# allocate coordinate array of the right size and type
# (the type float32 is crucial to match the underlying C-code!!)
x = np.zeros((natoms, DIM), dtype=np.float32)
# allocate unit cell box
box = np.zeros((DIM, DIM), dtype=np.float32)

# open file
XTC = xdrfile_open(xtc, 'r')

# loop through file until return status signifies end or a problem
# (it should become exdrENDOFFILE on the last iteration)
status = exdrOK
while status == exdrOK:
  status,step,time,prec = read_xtc(XTC, box, x)
  # do something with x
  centre = x.mean(axis=0)
  print 'Centre of geometry at %(time)g ps: %(centre)r' % vars()
  #print x

# finally close file
xdrfile_close(XTC)

