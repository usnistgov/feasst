import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as manimation
from matplotlib import pyplot
from shapely.geometry.polygon import LinearRing
from shapely.geometry.polygon import Polygon
from shapely.geometry.point import Point
from descartes import PolygonPatch
from numpy import linalg as LA
#from figures import SIZE
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar

#define shape parameters, units normalized such that the area is equivalent for any d
#for case emailed to me by John Royer, d=310nm, rg=60nm, L=1500nm, D=1661.251244, d/D=0.186606331, rg/D=0.036117354
rg=0.04
print("rg %.18f" % (rg))
#print 's%.18fs%.18fs%.18f' % "rg", rg, "d", d, "L", L

makeplot=1              # flag to make plot
makeplot=0              # flag to make plot
ngridd=51               # number of points for table in dist
ngridt=101              # number of points for table in angle
ngridr=2               # number of points for table in roundedness

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

def roundSquareOverlap(xd):
  """Return the overlap of two rounded squares"""
  # x: position of the center of second square
  x = np.arange(2, dtype=np.float)
  x[0] = xd
  x[1] = 0

  # create shapes
  x1 = [rs0, rs1, rs2, rs3]
  x2 = [rs02+x, rs12+x, rs22+x, rs32+x]
  ring = LinearRing(x1)
  ring2 = LinearRing(x2)
  sphere = Point(0, 0).buffer(ish)
  sphere2 = Point(x).buffer(ish)
  insquare = Polygon(x1)
  insquare2 = Polygon(x2)
  spheredep = Point(0, 0).buffer(ish+rg)
  spheredep2 = Point(x).buffer(ish+rg)
  roundsquare = ring.buffer(d, cap_style=1)
  roundsquare2 = ring2.buffer(d, cap_style=1)
  roundsquaredep = ring.buffer(d+rg, cap_style=1)
  roundsquaredep2 = ring2.buffer(d+rg, cap_style=1)
  dep = roundsquaredep.difference(roundsquare)
  dep2 = roundsquaredep2.difference(roundsquare2)
  tot = roundsquare.union(insquare)
  tot2 = roundsquare2.union(insquare2)
  totdep = roundsquaredep.union(spheredep)
  totdep2 = roundsquaredep2.union(spheredep2)
  overlap = tot2.intersection(tot)
  overlapdep = totdep2.intersection(totdep)
  da=overlap.area
  if (da > 0):
    da=1e115
  else:
    da=-overlapdep.area/vdep
  if (overlapdep.area == 0):
    da += xd*10
  #print "xd", xd, "da", da

  # plot
  if (makeplot == 1):
    ax = fig.add_subplot(111)
    #ax = fig.add_subplot(121)
    plot_coords(ax, ring)
    plot_line(ax, ring)
    plot_coords(ax, ring2)
    plot_line(ax, ring2)
    patch1 = PolygonPatch(tot, zorder=2)
    patch12 = PolygonPatch(tot2, zorder=2)
    #patch2 = PolygonPatch(roundsquaredep, zorder=2)
    patch2 = PolygonPatch(totdep, zorder=2)
    patch22 = PolygonPatch(totdep2, zorder=2)
    #patch1 = PolygonPatch(tot, zorder=2)
    #patch1 = PolygonPatch(roundsquare, zorder=2)
    #patch1 = PolygonPatch(roundsquare, alpha=0.5, zorder=2)
    ax.add_patch(patch2)
    ax.add_patch(patch22)
    ax.add_patch(patch1)
    ax.add_patch(patch12)
    xrange = [-1, 2]
    yrange = [-1, 1]
    ax.set_xlim(*xrange)
    ax.set_xticks(range(*xrange) + [xrange[-1]])
    ax.set_ylim(*yrange)
    ax.set_yticks(range(*yrange) + [yrange[-1]])
    ax.set_aspect(1)
    pyplot.savefig('cubed'+str(d)+'rg'+str(rg)+'.pdf')
    pyplot.savefig('cubed'+str(d)+'rg'+str(rg)+'.jpg')
    writer.grab_frame()
    pyplot.clf()
  return da

def roundSquareOverlapSq(xd):
  return roundSquareOverlap(xd)**2

fig = pyplot.figure(1, dpi=9000, figsize=(5,5))

#compute maximum and minimium distances for all range of d [0, 0.5]
absoluteMinDist = math.sqrt(math.pi/4)
absoluteMaxDist = 2*rg+math.sqrt(math.pi/2)
ngriddt=(math.pi/2)/float(ngridt-1)
ngridth=int(ngridt/2)+1
ngridrmax=0.5
print ngridth
app="rg"+str(rg)+"nt"+str(ngridt)+"nz"+str(ngridd)+"nd"+str(ngridr)
rm=open("rm"+app, 'w')
rc=open("rc"+app, 'w')
vt=open("vt"+app, 'w')
rst=open("rst"+app, 'w')

rst.write("# className PairTabular\n")
rst.write("# rCut %.16f\n" % absoluteMaxDist)
rst.write("# rCutInner %.16f\n" % absoluteMinDist)
rst.write("# tabHardFileName /home/hwh/hwhcms/table/roundSquare/rm"+app+"\n")
rst.write("# tabCutFileName /home/hwh/hwhcms/table/roundSquare/rc"+app+"\n")
rst.write("# tabPEFileName /home/hwh/hwhcms/table/roundSquare/vt"+app+"\n")

rm.write("# tabType roundSquareMin\n")
rm.write("# rg %.16f\n" % rg)
rm.write("# rCut %.16f\n" % absoluteMaxDist)
rm.write("# rCutInner %.16f\n" % absoluteMinDist)
rm.write("# tabDims 3\n")
rm.write("# dimorder0 2\n")
rm.write("# dimn0 %i\n" % ngridr)
rm.write("# dimmin0 0\n")
rm.write("# dimmax0 %.16f\n" % ngridrmax)
rm.write("# dimorder1 0\n")
rm.write("# dimn1 %i\n" % ngridth)
rm.write("# dimmin1 0\n")
rm.write("# dimmax1 %.16f\n" % float(math.pi/4.))
rm.write("# dimorder2 1\n")
rm.write("# dimn2 %i\n" % ngridt)
rm.write("# dimmin2 %.16f\n" % float(-math.pi/4))
rm.write("# dimmax2 %.16f\n" % float(math.pi/4))

rc.write("# tabType roundSquareCut\n")
rc.write("# rg %.16f\n" % rg)
rc.write("# rCut %.16f\n" % absoluteMaxDist)
rc.write("# rCutInner %.16f\n" % absoluteMinDist)
rc.write("# tabDims 3\n")
rc.write("# dimorder0 2\n")
rc.write("# dimn0 %i\n" % ngridr)
rc.write("# dimmin0 0\n")
rc.write("# dimmax0 %.16f\n" % ngridrmax)
rc.write("# dimorder1 0\n")
rc.write("# dimn1 %i\n" % ngridth)
rc.write("# dimmin1 0\n")
rc.write("# dimmax1 %.16f\n" % float(math.pi/4.))
rc.write("# dimorder2 1\n")
rc.write("# dimn2 %i\n" % ngridt)
rc.write("# dimmin2 %.16f\n" % float(-math.pi/4))
rc.write("# dimmax2 %.16f\n" % float(math.pi/4))

vt.write("# tabType roundSquareDep\n")
vt.write("# rg %.16f\n" % rg)
vt.write("# rCut %.16f\n" % absoluteMaxDist)
vt.write("# rCutInner %.16f\n" % absoluteMinDist)
vt.write("# tabDims 4\n")
vt.write("# dimorder0 3\n")
vt.write("# dimn0 %i\n" % ngridr)
vt.write("# dimmin0 0\n")
vt.write("# dimmax0 %.16f\n" % ngridrmax)
vt.write("# dimorder1 0\n")
vt.write("# dimn1 %i\n" % ngridd)
vt.write("# dimmin1 0\n")
vt.write("# dimmax1 1\n")
vt.write("# dimorder2 1\n")
vt.write("# dimn2 %i\n" % ngridth)
vt.write("# dimmin2 0\n")
vt.write("# dimmax2 %.16f\n" % float(math.pi/4))
vt.write("# dimorder3 2\n")
vt.write("# dimn3 %i\n" % ngridt)
vt.write("# dimmin3 %.16f\n" % float(-math.pi/4))
vt.write("# dimmax3 %.16f\n" % float(math.pi/4))

# define rotation matrix for pi/2 rotational symmetry of square
rotmat = np.arange(4, dtype=np.float).reshape((2,2))
rotmat[0,0] = 0
rotmat[0,1] = -1
rotmat[1,0] = 1
rotmat[1,1] = 0

#writer.setup(fig, 'cubed'+str(d)+'rg'+str(rg)+'.mp4', dpi=90)
with writer.saving(fig, 'cubed'+'rg'+str(rg)+'.mp4', ngridd*ngridt*ngridth*ngridr):
  # r1: orientation vector, which is normalized to point in the direction of a corner
  r = np.arange(2, dtype=np.float)
  unit = np.arange(2, dtype=np.float)
  unit[0] = 1
  unit[1] = 0

  # r2: orientation vector, which is normalized to point in the direction of a corner
  r2 = np.arange(2, dtype=np.float)
  #for rd2 in np.arange(0, 0.5, 1):
  #for theta2 in np.arange(0, (math.pi/4)+ngriddt*0.5, ngriddt):
  dngriddr=ngridrmax/(ngridr-1)
  for d in np.arange(0, ngridrmax+0.5*dngriddr, dngriddr):
    # compute derived quantities and print
    L=math.sqrt( (math.pi/4)+(4-math.pi)*d**2)
    vdep=math.pi*rg**2
    innersquare=L-2*d       # side length of square inscribed by circle
    ish=0.5*innersquare     # half of the side length of square inscribed by circle
    maxdist=2*(d+rg+ish*math.sqrt(2))   #largest cutoff distance for corner-corner interaction
  
    for theta2 in np.arange(-math.pi/4, (math.pi/4)+ngriddt*0.5, ngriddt):
      r2[0] = math.cos(theta2)
      # note, for second angle, defined theta only rij axis, so y value is negative
      r2[1] = -math.sin(theta2)
      rs02=r2*ish*math.sqrt(2)
      rs12=np.dot(rotmat, rs02)
      rs22=np.dot(rotmat, rs12)
      rs32=np.dot(rotmat, rs22)

      #for theta in np.arange(0, 0.5, 1):
      for theta in np.arange(0, (math.pi/4)+ngriddt*0.5, ngriddt):
      #for theta in np.arange(-math.pi/4, (math.pi/4)+ngriddt*0.5, ngriddt):
        r[0] = math.cos(theta)
        r[1] = math.sin(theta)
        rs0=r*ish*math.sqrt(2)
        rs1=np.dot(rotmat, rs0)
        rs2=np.dot(rotmat, rs1)
        rs3=np.dot(rotmat, rs2)
        eps=0.00001
        resrm = minimize_scalar(roundSquareOverlap, method="Brent", bracket=(L-eps, maxdist, maxdist+eps))
        #print resrm.x, maxdist, maxdist+eps
        #print roundSquareOverlapSq(resrm.x), roundSquareOverlapSq(maxdist), roundSquareOverlapSq(maxdist+eps)
        resrc = minimize_scalar(roundSquareOverlapSq, method="Brent", bracket=(resrm.x - eps, maxdist, maxdist+eps))
        print resrm.fun, "dx", resrm.x, "dcut", resrc.x, "t1", theta, "t2", theta2
        rm.write("%.8f d %f t1 %.16f t2 %.16f\n" % (resrm.x, d, theta, theta2))
        #print >>rm, resrm.x, "t1", theta, "t2", theta2
        rc.write("%.8f d %f t1 %.16f t2 %.16f\n" % (resrc.x, d, theta, theta2))
        #print >>rc, resrc.x, "t1", theta, "t2", theta2

        #now generate table between rm and rc
        ngriddd=1/float(ngridd-1)
        #ngriddd=(resrc.x - resrm.x)/float(ngridd-1)
        for z in np.arange(0, 1+0.5*ngriddd, ngriddd):
          vt.write("%.8f d %f z %f t1 %.16f t2 %.16f\n" % (roundSquareOverlap(z*(resrc.x - resrm.x) + resrm.x), d, z, theta, theta2))
          #print >>vt, roundSquareOverlap(z*(resrc.x - resrm.x) + resrm.x), "z", z, "t1", theta, "t2", theta2

