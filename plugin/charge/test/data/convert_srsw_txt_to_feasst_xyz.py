"""
This script was used to convert from SRSW txt file configuration to FEASST-readable xyz file.
"""

import numpy as np
import pandas as pd
from pyfeasst import physical_constants

for U in [931.15451, -34.16569, 371.46525, -6046.43627, 95078.89447, -96297.75579]: # kJ/mol
    print('U(kJ/mol)', U, 'U(K)', U*1e3/physical_constants.MolarGasConstant().value())

for U in [2.50251E+04, -1.63091E+02, 2.23372E+04, -1.71462E+05, 2.85884E+06, -2.89549E+06]: # K
    print('U(K)', U, 'U(kJ/mol)', U*physical_constants.MolarGasConstant().value()/1e3)

print('U_ewald_config4', 185.60955025935957*1e3/physical_constants.MolarGasConstant().value())
print('U_ewald_config4', 185.720351510974*1e3/physical_constants.MolarGasConstant().value())

lx=30
ly=28.97777478867205
lz=29.51512917398008
xy=7.764571353075622
xz=-2.6146722824297473
yz=-4.692615336756641

a=lx
b=np.sqrt(ly*ly+xy*xy)
c=np.sqrt(lz*lz+xz*xz+yz*yz)
print('abc', a, b, c)
print('alpha', np.arccos((xy*xz+ly*yz)/b/c)*180/np.pi)
print('beta', np.arccos(xz/c)*180/np.pi)
print('gamma', np.arccos(xy/b)*180/np.pi)

A=30
B=30
C=30
alpha=100/180*np.pi
beta=95/180*np.pi
gamma=75/180*np.pi

A=36
B=36
C=36
alpha=90/180*np.pi
beta=60/180*np.pi
gamma=90/180*np.pi

lx = A
xy = B*np.cos(gamma)
xz = C*np.cos(beta)
ly=np.sqrt(B*B-xy*xy)
yz=(B*C*np.cos(alpha)-xy*xz)/ly
lz=np.sqrt(C*C-xz*xz-yz*yz)
print('l', lx, ly, lz)
print('xy, xz, yz', xy, xz, yz)

#df = pd.read_csv('spce_triclinic_sample_periodic1.xyz', header=None, skiprows=2, delim_whitespace=True)
df = pd.read_csv('spce_monoclinic_sample_periodic4.xyz', header=None, skiprows=2, delim_whitespace=True)
print(df)

num_spce = int(len(df)/3)
print('num_spce', num_spce)
#for line in len(df):
theta_HOH = 109.47/180*np.pi
bond_OH = 1.
lx=36
ly=36
lz=31.17691453623979
xy=0.
xz=18.
yz=0.

for spce in range(num_spce):
    ox = df[1][spce*3]
    oy = df[2][spce*3]
    oz = df[3][spce*3]
    h1x = df[1][spce*3+1]
    h1y = df[2][spce*3+1]
    h1z = df[3][spce*3+1]
    h2x = df[1][spce*3+2]
    h2y = df[2][spce*3+2]
    h2z = df[3][spce*3+2]
    dh1x=h1x-ox
    dh1y=h1y-oy
    dh1z=h1z-oz
    dh2x=h2x-ox
    dh2y=h2y-oy
    dh2z=h2z-oz
    h1_bond_length = np.sqrt(dh1x*dh1x+dh1y*dh1y+dh1z*dh1z)
    h2_bond_length = np.sqrt(dh2x*dh2x+dh2y*dh2y+dh2z*dh2z)
    if np.abs(h1_bond_length - 1) > 1e-8:
        #print(spce*3+4, 'h1_bond_length', h1_bond_length, dh1x, dh1y, dh1z)
        if dh1z > 1.1*bond_OH:
            dh1z -= lz
            dh1y -= yz
            dh1x -= xz
        elif dh1z < -1.1*bond_OH:
            dh1z += lz
            dh1y += yz
            dh1x += xz
        elif dh1y > 1.1*bond_OH:
            assert False
        elif dh1x > 1.1*bond_OH:
            dh1x -= lx
        elif dh1x < -1.1*bond_OH:
            dh1x += lx
        h1_bond_length = np.sqrt(dh1x*dh1x+dh1y*dh1y+dh1z*dh1z)
        #print(spce*3+4, 'fixed h1_bond_length', h1_bond_length, dh1x, dh1y, dh1z)
        print('line', spce*3+4, 'xyz', ox+dh1x, oy+dh1y, oz+dh1z)
    if np.abs(h2_bond_length - 1) > 1e-8:
        #print(spce*3+5, 'h2_bond_length', h2_bond_length, dh2x, dh2y, dh2z)
        if dh2z > 1.1*bond_OH:
            dh2z -= lz
            dh2y -= yz
            dh2x -= xz
        elif dh2z < -1.1*bond_OH:
            dh2z += lz
            dh2y += yz
            dh2x += xz
        elif dh2y > 1.1*bond_OH:
            assert False
        elif dh2x > 1.1*bond_OH:
            dh2x -= lx
        elif dh2x < -1.1*bond_OH:
            dh2x += lx
        h2_bond_length = np.sqrt(dh2x*dh2x+dh2y*dh2y+dh2z*dh2z)
        #print(spce*3+5, 'fixed h2_bond_length', h2_bond_length, dh2x, dh2y, dh2z)
        print('line', spce*3+5, 'xyz', ox+dh2x, oy+dh2y, oz+dh2z)

