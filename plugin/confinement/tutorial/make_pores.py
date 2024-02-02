import os
import copy
import argparse
import numpy as np
basename = os.path.basename(__file__)[:-3]

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--spacing', type=float, default=9./8., help='spacing between each site')
PARSER.add_argument('--size', type=int, default=4., help='number of sites in a pore in each dimension')
PARSER.add_argument('--show_plot', type=int, default=1., help='if 1, show plot of one pore')
# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['file_name'] = basename + 'd' + str(PARAMS['spacing']) + 's' + str(PARAMS['size']) + '.fstprt'

xs = list()
ys = list()
zs = list()
for x in np.arange(0.5, PARAMS['size'], 1):
  for y in np.arange(0.5, PARAMS['size'], 1):
    for z in np.arange(0.5, PARAMS['size'], 1):
      if (np.abs(x) < 1 or np.abs(x) > PARAMS['size']-1 or \
          np.abs(y) < 1 or np.abs(y) > PARAMS['size']-1 or \
          np.abs(z) < 1 or np.abs(z) > PARAMS['size']-1) and \
         (x != PARAMS['size']-1.5 or y != PARAMS['size']-1.5) and \
         (x != PARAMS['size']-1.5 or z != PARAMS['size']-1.5) and \
         (y != PARAMS['size']-1.5 or z != PARAMS['size']-1.5):
        xs.append(x)
        ys.append(y)
        zs.append(z)

xs = np.array(xs)
ys = np.array(ys)
zs = np.array(zs)
xs *= PARAMS['spacing']
ys *= PARAMS['spacing']
zs *= PARAMS['spacing']
xsa=np.zeros(0)
ysa=np.zeros(0)
zsa=np.zeros(0)
for xshift in [0, PARAMS['size']*PARAMS['spacing']]:
  for yshift in [0, PARAMS['size']*PARAMS['spacing']]:
    for zshift in [0, PARAMS['size']*PARAMS['spacing']]:
      xsa = np.concatenate((xsa, xs + xshift))
      ysa = np.concatenate((ysa, ys + yshift))
      zsa = np.concatenate((zsa, zs + zshift))

f = open(PARAMS['file_name'], 'w')
f.write('# 2x2x2 pores\n\
\n\
'+str(len(xsa))+' sites\n\
\n\
1 site types\n\
\n\
Site Properties\n\
\n\
0 sigma 1.0 epsilon 1.0 cutoff 2.5\n\
\n\
Sites\n\
\n\
')

for ix,x in enumerate(xsa):
  f.write('0 0 '+str(x) + ' ' + str(ysa[ix]) + ' ' + str(zsa[ix]) + '\n')

if PARAMS['show_plot']:
  import matplotlib.pyplot as plt
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
  ax.scatter(xs, ys, zs, s=10)
  #ax.scatter(xsa, ysa, zsa, s=10)
  plt.show()
