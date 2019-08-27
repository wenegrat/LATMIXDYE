#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate maximum displacements

Created on Mon Mar 26 09:44:57 2018

@author: jacob
"""

import numpy as np
import h5py
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import scipy.stats as stats
from astropy.stats import circstats
#%% LOAD
filename = '/home/jacob/dedalus/LATMIX/front.mat' # Pick which run
f = h5py.File(filename, 'r')

# List all groups
print("Keys: %s" % f.keys())
keys = list(f.keys())[:]

ua = f['ua']



dye1 = f['dye1']
x = f['x'][:,0]
z = f['z'][:,0]
z = z[1:]
time = f['t'][:,0]
dye1 = dye1[:,1:,:]
nt, nz, nx = dye1.shape

#%% Limit analysis to period where dye is released
indrelease = 162
indend = np.argmin(np.abs(time/86400 + 64.5 - 66))

zu = np.argmin(np.abs(z+5))
zd100 = np.argmin(np.abs(z+100))
zd50 = np.argmin(np.abs(z+50))
#%% First calculate maximum displacements by averaged u velocity
mu100 = integrate.trapz(ua[:,zd100:zu], x = z[zd100:zu], axis=-1)/(z[zu]-z[zd100]) # might try restricting this to the upper 100 meter or so
up100 = ua - mu100[:,np.newaxis]
mu50 = integrate.trapz(ua[:,zd50:zu], x = z[zd50:zu], axis=-1)/(z[zu]-z[zd50]) # might try restricting this to the upper 100 meter or so
up50 = ua - mu50[:,np.newaxis]
xdis100 = np.zeros((time.size, nz))
xdis50 = np.zeros((time.size, nz))
for i in range(zd100, zu):
    xdis100[:,i] = integrate.cumtrapz(up100[0:1569,i], x= time, initial=0)

for i in range(zd50, zu):
    xdis50[:,i] = integrate.cumtrapz(up50[0:1569,i], x= time, initial=0)
#%% Calculate X displacements at each depth relative to indrelease
xdis100 = xdis100 - xdis100[indrelease,:]
xdis100[0:indrelease-1,:] = 0

maxdelta100 = np.max(xdis100, axis=-1) - np.min(xdis100, axis=-1)

xdis50 = xdis50 - xdis50[indrelease,:]
xdis50[0:indrelease-1,:] = 0

maxdelta50 = np.max(xdis50, axis=-1) - np.min(xdis50, axis=-1)
#%%
plt.figure()
plt.plot(time/86400 + 64.5, maxdelta100, label='z=(-5, -100)')
plt.plot(time/86400 + 64.5, maxdelta50 , label='z=(-5, -50)')
plt.xlabel('yearday')
plt.legend()
plt.ylabel('Maximum $\Delta x$')
plt.tight_layout()
plt.grid()
#%%
plt.figure()
plt.pcolor(time/86400, z, np.transpose(xdis))
plt.colorbar()
