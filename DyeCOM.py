#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
COM Statistics

Created on Mon Mar 26 09:44:57 2018

@author: jacob
"""

import numpy as np
import h5py
import scipy.integrate as integrate
import matplotlib.pyplot as plt
#%% LOAD
filename = '/home/jacob/dedalus/LATMIX/front.mat' # Pick which run
f = h5py.File(filename, 'r')

# List all groups
print("Keys: %s" % f.keys())
keys = list(f.keys())[:]

u = f['u']
w = f['w']

dye1 = f['dye1']
dye2 = f['dye2']
time = f['t'][:,0]
b = f['b'] # total buoyancy
bp = f['bp'] # perturbation buoyancy
x = f['x'][:,0]
z = f['z'][:,0]

nt, nz, nx = u.shape
#%% COM Calcs
ubar = integrate.trapz(u[:,0:-1,:],x=z[0:-1], axis=1)/z[0]

# define circular coordinate system for averaging
theta = x/np.max(x)*2*np.pi
xi = np.cos(theta)
zi = np.sin(theta)

# Pick which variable to use
dyevar = dye1

# Pre-allocate
zeta = np.zeros((nt, nx))
TT = np.zeros((nt,))
ZCM = np.zeros((nt,))
xposcom = np.zeros((nt,))
zetas = np.zeros((nt, nx))
xs = np.zeros((nt, nx))
ubars = np.zeros((nt, nx))
XVar = np.zeros((nt,))

for i in range(0, nt):
    
    zeta[i,:] = integrate.trapz(dyevar[i,1:,:], x=z[1:], axis=0) # vertically integrated dye content
    TT[i] = integrate.trapz(zeta[i,:], x=x) # Total domain tracer content (shouldn't this remain constant?)
    
    # Calculate zonal center of mass
    xii = xi*zeta[i,:]
    zii = zi*zeta[i,:]
    tbar = np.arctan2(-np.mean(zii), -np.mean(xii)) + np.pi
    ZCM[i] = np.max(x)*tbar/(2*np.pi) 
    
    # Find index of Zonal Center of Mass
    ind = np.where(x>=ZCM[i])[0][0]
    xposcom[i] = x[ind]
    indshift = ind-nx/2 # This is used to center the dye in the domain.
    indshift = indshift.astype('int64')
    
    #Circular shfit everything
    zetas[i,:] = np.roll(zeta[i,:], -indshift, axis=0)
    #xs[i,:] = np.roll(x, -indshift)
    xs = x - x[int(nx/2)]
    ubars[i,:] = np.roll(ubar[i,:], -indshift, axis=0)
    
    # Calculate dye variance around COM
    M20 = integrate.trapz(xs**2*zetas[i,:], x=x, axis=0) # check this, should be xs? should x = 0 at COM?
    XVar[i] = M20/TT[i] # maybe this sets x=0 at center? Not important for Kappa Calc....
    
xposunwrap = np.unwrap(xposcom*2*np.pi/x[-1] - np.pi)*x[-1]/(2*np.pi)
#%%
# Calculate kappa and buoyancy class
kappa = 0.5*np.gradient(XVar[0:time.size])/np.gradient(time)

#%%
plt.figure()
plt.plot(time/86400+64.5,kappa[0:time.size])
plt.xlabel('Yearday')
plt.ylabel('$\kappa_h$ m$^2$s$^{-1}$')
plt.grid()
plt.tight_layout()
plt.xlim((64.5, 66))

#%% VARIANCE
plt.figure()
plt.plot(time/86400+64.5,XVar[0:time.size])
plt.xlabel('Yearday')
plt.ylabel('$\kappa_h$ m$^2$s$^{-1}$')
plt.grid()
plt.tight_layout()
plt.xlim((64.5, 66))

#%%
plt.figure()
plt.plot(xs, np.transpose(zetas[200:900:50,:]))
