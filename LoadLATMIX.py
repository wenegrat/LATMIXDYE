#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 15:20:39 2018

@author: jacob
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
#mpl.verbose.set_level("helpful")
#import matplotlib.pyplot as plt

#from matplotlib.gridspec import GridSpec
#matplotlib.use("Agg")
mpl.pyplot.rcParams['animation.mencoder_path'] = '/usr/bin/mencoder'

import matplotlib.animation as animation
import h5py
import scipy.io as spio

#%% LOAD
filename = '/home/jacob/dedalus/LATMIX/front.mat' # Pick which run
f = h5py.File(filename, 'r')

# List all groups
print("Keys: %s" % f.keys())
keys = list(f.keys())[:]

u = f['u']
w = f['w']
q = f['q']
qf = f['qv']
qv = f['q_z']
qh = f['q_h']
dye1 = f['dye1']
dye2 = f['dye2']
time = f['t'][:,0]
b = f['b'] # total buoyancy
bp = f['bp'] # perturbation buoyancy
x = f['x'][:,0]
z = f['z'][:,0]

uflux1 = f['uflux1']
wflux1 = f['wflux1']
nt, nz, nx = u.shape

#%% PV Plot
tn = 1100
fig, ax = plt.subplots(3,1)
IM = ax[0].pcolor(time/86400, z[0:-2], np.transpose(q[0, :,0:-2]), vmin = -1e-8, vmax=1e-8)
plt.colorbar(IM, ax=ax[0])

IM = ax[1].pcolor(time/86400, z[0:-2], np.transpose(qv[0, :,0:-2]), vmin = -1e-8, vmax=1e-8)
plt.colorbar(IM, ax=ax[1])

IM = ax[2].pcolor(time/86400, z[0:-2], np.transpose(q[0, :,0:-2]-qv[0, :,0:-2]), vmin = -1e-8, vmax=1e-8)
plt.colorbar(IM, ax=ax[2])
#%%
var = dye1
bc = np.linspace(0, 0.01, 20)
fig = mpl.pyplot.figure(figsize=(12, 5))
ax1 = mpl.pyplot.axes()

bc = np.linspace(0, 0.01, 20)

#plt.hold(True)
#We need to prime the pump, so to speak and create a quadmesh for plt to work with
quad1 = ax1.pcolormesh(x, z, var[0,:,:])
quad2 = ax1.contour(x, z, b[0,:,:], bc)

#quad1.set_array(dye1[100,:,:].ravel())
def animate(i):
    global quad2
    tn = time[i]/86400
    mpl.pyplot.title('Time: %.2f'%tn)
    print(i)
    dyn = var[i,:,:]/np.max(var[i,:,:])
    #This is where new data is inserted into the plot.
    quad1.set_array(dyn.ravel())
#    quad2.set_array(b[i,:,:].ravel())
    for c in quad2.collections:
        c.remove()
    quad2 = ax1.contour(x, z, b[i,:,:], bc)

    return quad2, quad1

anim = mpl.animation.FuncAnimation(fig, animate, frames = 1400, blit = False)