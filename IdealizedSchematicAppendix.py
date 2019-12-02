#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 13:35:03 2019

@author: jacob
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cmocean.cm as cmo
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.axis3d import Axis
if not hasattr(Axis, "_get_coord_info_old"):
    def _get_coord_info_new(self, renderer):
        mins, maxs, centers, deltas, tc, highs = self._get_coord_info_old(renderer)
        mins += deltas / 4
        maxs -= deltas / 4
        return mins, maxs, centers, deltas, tc, highs
    Axis._get_coord_info_old = Axis._get_coord_info  
    Axis._get_coord_info = _get_coord_info_new
plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 32
plt.rcParams['contour.negative_linestyle'] = 'solid'
#%%

z = np.linspace(-100, 0, 100)

UZI = (z-z[0])*0.1
VZI = 0*z

UZM = UZI.copy()
zind = np.argmin(np.abs(z+50))

UZM[zind:] = np.mean(UZI[zind:])

time = np.linspace(0, 2*np.pi, 100)

nt = time.size
nz = z.size
u = np.zeros((nt, nz))
v = np.zeros((nt, nz))


ueq = UZM - UZI #ueq is difference between geo solution and final state
for i in range(0, nz):
    u[:,i] = UZI[i] + ueq[i]*(1-np.cos(time))
    v[:,i] = ueq[i]*np.sin(time)
#%%
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

ax.plot(UZI, VZI, z, color='k', linestyle='solid', linewidth=2)
ax.plot(UZM, VZI, z, color='k', linewidth=2, linestyle='dashed')


for j in range(zind, nz, 5):
#    ax.plot(u[:,j], v[:,j], z[j])
    for i in range(nt-1):
        ax.plot(u[i:i+2, j], v[i:i+2,j], z[j], color=cmo.phase(i/nt))
    
al = 10
ax.set_xlim((0, al))
ax.set_ylim((-al/2, al/2))
ax.set_xticklabels([0])
ax.set_yticklabels(['','',0,'',''], rotation = 0, va='center', ha='left')
ax.set_zticklabels([])
ax.set_xlabel('u')
ax.set_ylabel('v', labelpad=10)
ax.set_zlabel('z')
ax.view_init(elev = 25, azim =-65)

sm = plt.cm.ScalarMappable(cmap=cmo.phase, norm=plt.Normalize(vmin=0, vmax=2*np.pi))
# fake up the array of the scalar mappable. Urgh...
sm._A = []

#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(sm, fraction=0.03, pad=0.03)
cb.set_label('Normalized time (ft)')
cb.set_ticks([0, 2*np.pi])
cb.set_ticklabels(['0', '$2\pi$'])
ax.set_zlim((z[0], z[-1] ))
plt.tight_layout()

#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/FigureS2.pdf', bbox_inches='tight')