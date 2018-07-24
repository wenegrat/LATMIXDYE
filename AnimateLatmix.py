#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AnimateLatmix.py
Created on Mon Mar 26 09:43:29 2018

@author: jacob
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
mpl.pyplot.rcParams['animation.mencoder_path'] = '/usr/bin/mencoder'

import matplotlib.animation as animation
import h5py
import scipy.io as spio
#%% LOAD NL SIM SNAPS
filename = '/home/jacob/dedalus/LATMIX/front.mat'
f = h5py.File(filename, 'r')

# List all groups
print("Keys: %s" % f.keys())
keys = list(f.keys())[:]

u = f['u']
q = f['q'][0,:,:]
#v = f['v']
w = f['w']
tx = f['tx']
ty = f['ty']
dye1 = f['dye1']
dye2 = f['dye2']
time = f['t'][:,0]
b = f['b'] # total buoyancy
bp = f['bp'] # perturbation buoyancy
x = f['x'][:,0]
z = f['z'][:,0]

#%% SAVE MOVIE

Writer = mpl.animation.writers['mencoder']
writer = Writer(fps=10, metadata=dict(artist='Me'),extra_args=["-really-quiet"], bitrate=2000)

var = u
bc = np.linspace(0, 0.01, 20)
fig = mpl.pyplot.figure(figsize=(12, 5))
ax1 = mpl.pyplot.axes()

bc = np.linspace(0, 0.01, 20)

#We need to prime the pump, so to speak and create a quadmesh for plt to work with
quad1 = ax1.pcolormesh(x, z, var[0,:,:], vmin=-0.2, vmax=0.2)
quad2 = ax1.contour(x, z, b[0,:,:], bc)

def animate(i):
    global quad2
    tn = time[i]/86400
    mpl.pyplot.title('Time: %.2f'%tn)
    print(i)
    dyn = var[i,:,:]#/np.max(var[i,:,:])
    #This is where new data is inserted into the plot.
    quad1.set_array(dyn.ravel())
#    quad2.set_array(b[i,:,:].ravel())
    for c in quad2.collections:
        c.remove()
    quad2 = ax1.contour(x, z, b[i,:,:], bc)

    return quad2, quad1

# Frames argument is how many time steps to call 'animate' function
anim = mpl.animation.FuncAnimation(fig, animate, frames = 450, blit = False)
#anim.save('u_anim.avi', writer=writer) # Uncomment to Save Animation

#mpl.pyplot.show() # Can turn this on to just show the animation on screen
print('finished')


