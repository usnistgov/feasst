import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as manimation
from matplotlib import pyplot
from shapely.geometry.polygon import LinearRing
from shapely.geometry.point import Point
from descartes import PolygonPatch
from numpy import linalg as LA

# define shape parameters, units normalized by side length
filename = 'movie.xyz'  # xyz movie file name
d=(310/float(750)/2)    # radius of rounded corners
rg=60/float(750)/2      # radius of depletant
skipevery=100           # skip frames

# compute derived quantities and print
innersquare=1-2*d       # side length of square inscribed by circle
ish=0.5*innersquare     # half of the side length of square inscribed by circle
maxdist=2*(d+rg+ish*math.sqrt(2))

#precompute
rotmat = np.arange(4, dtype=np.float).reshape((2,2))
rotmat[0,0] = 0
rotmat[0,1] = -1
rotmat[1,0] = 1
rotmat[1,1] = 0
rij = np.arange(2, dtype=np.float)
x = np.arange(2, dtype=np.float)

fig = pyplot.figure(1, dpi=500, figsize=(600,600))
#fig = pyplot.figure(1, dpi=90, figsize=(100,100))
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
        comment='Movie support!')
writer = FFMpegWriter(fps=15, metadata=metadata)

COLOR = {
    True:  '#6699cc',
    False: '#ff3333'
    }

def v_color(ob):
    return COLOR[ob.is_valid]

def plot_coords(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, 'o', color='#999999', zorder=1)

def plot_line(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, color=v_color(ob), alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)

#determine number of frames and lines
f = open(filename)
lines = f.readlines()
f.close()
linenum=0
num_lines = sum(1 for line in open(filename))
nframes = 0
iframe = 0
while (linenum < num_lines):
  nMol=int(lines[linenum])/2
  linenum = linenum + 2 + 2*nMol
  if (iframe % skipevery == 0):
    nframes = nframes + 1
  iframe += 1
linenum=0
print "nframes", nframes*skipevery


iframe = 0
#writer.setup(fig, 'cubed'+str(d)+'rg'+str(rg)+'.mp4', dpi=90)
#with writer.saving(fig, 'cubed'+str(d)+'rg'+str(rg)+'.mp4',nframes+1):
with writer.saving(fig, 'movied'+str(round(d,4))+'rg'+str(rg)+'.mp4', 3):
  while (linenum < num_lines):
    nMol=int(lines[linenum])/2
    linenum = linenum + 2
    #print iframe % skipevery
    if (iframe % skipevery != 0):
      linenum += nMol*2
      iframe += 1
    else:
      ax = fig.add_subplot(111)
      xrange = [-6, 6]
      yrange = [-6, 6]
      ax.set_xlim(*xrange)
      ax.set_xticks(range(*xrange) + [xrange[-1]])
      ax.set_ylim(*yrange)
      ax.set_yticks(range(*yrange) + [yrange[-1]])
      ax.set_aspect(1)
      for iatom in range(nMol):
        columns = lines[linenum].split(' ')
        x[0] = float(columns[1])
        x[1] = float(columns[2])
        #print x
        linenum = linenum + 1
        columns = lines[linenum].split(' ')
        x2 = float(columns[1])
        y2 = float(columns[2])
        #print x2, y2
        linenum = linenum + 1
        rij[0] = x[0] - x2
        rij[1] = x[1] - y2
        #print rij
        rs0=rij*ish*math.sqrt(2)
        rs1=np.dot(rotmat, rs0)
        rs2=np.dot(rotmat, rs1)
        rs3=np.dot(rotmat, rs2)
        ring = LinearRing([rs0+x, rs1+x, rs2+x, rs3+x])
        sphere = Point(x).buffer(ish)
        spheredep = Point(x).buffer(ish+rg)
        roundsquare = ring.buffer(d, cap_style=1)
        roundsquaredep = ring.buffer(d+rg, cap_style=1)
        dep = roundsquaredep.difference(roundsquare)
        tot = roundsquare.union(sphere)
        totdep = roundsquaredep.union(spheredep)
        plot_coords(ax, ring)
        plot_line(ax, ring)
        patch1 = PolygonPatch(tot, zorder=2)
        patch2 = PolygonPatch(totdep, zorder=2)
        ax.add_patch(patch2)
        ax.add_patch(patch1)
        #ax.set_title('a) valid')
      #pyplot.savefig('cubed'+str(d)+'rg'+str(rg)+'.pdf')
      #pyplot.savefig('cubed'+str(d)+'rg'+str(rg)+'.jpg')
      writer.grab_frame()
      pyplot.clf()
      iframe = iframe + 1
      print round(iframe/float(nframes),3), "percent complete: iframe", iframe, "of nframes", nframes*skipevery


