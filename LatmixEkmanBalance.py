#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 14:51:40 2018

@author: jacob
"""
import h5py

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
import scipy.integrate as integrate
import scipy.interpolate as interpolate
#%% LOAD NL SIM SNAPS
filename = '/home/jacob/dedalus/LATMIX/front.mat'
f = h5py.File(filename, 'r')

# List all groups
print("Keys: %s" % f.keys())
keys = list(f.keys())[:]

vw = f['vw']
uw = f['uw']
u = f['u']
#v = f['v']
vm = f['va']
tx = f['tx']
ty = f['ty']
#z = f['z'][:,0]
time = f['t'][:,0]
um = np.mean(u, axis=-1)
vwz = np.gradient(vw, axis=-1)/np.gradient(z)
uwz = np.gradient(uw, axis=-1)/np.gradient(z)

coriolis = 9.3e-5
#ue = -vwz/coriolis




#%% SETUP TIME DEPENDENT PROBLEM

deltat = np.diff(time)
ue = np.zeros((deltat.size, z.size))
ve = np.zeros((deltat.size, z.size))

ue[0,:] = um[0,:]
ve[0,:] = vm[0,:] 
for i in range(1, deltat.size):
    ue[i,:] = ue[i-1,:] + (coriolis*ve[i-1,:] - uwz[i,:])*deltat[i]
    if i<1000:
        print(i)
        ve[i,:] = ve[i-1,:] + (-coriolis*ue[i-1,:] - vwz[i,:])*deltat[i]
    else:
        ve[i,:] = ve[i-1,:] + (-coriolis*ue[i-1,:] - 0*vwz[i,:])*deltat[i]
#%%
tn = np.where((time/86400>0.5))[0][0]
plt.figure()
plt.plot(ue[tn,0:-2], z[0:-2])
plt.plot(um[tn,0:-2], z[0:-2])
plt.title(time[tn]/86400)
#plt.plot(np.mean(um[0:tn,0:-2], axis=0), z[0:-2])
    
#%%
fig, ax = plt.subplots(2, 1, sharex=False)
tl = 200

vmin = -0.15
vmax = -vmin
ax[0].pcolor(time[0:tl]/86400, z[0:-2], np.transpose(um[0:tl,0:-2]), vmin=vmin, vmax = vmax)
ax[0].set_title('Modeled u')
ax[1].pcolor(time[0:tl]/86400, z[0:-2], np.transpose(ue[0:tl,0:-2]), vmin=vmin, vmax = vmax)
ax[1].set_title('Ekman u (v stress only)')
plt.tight_layout()
ax[0].set_ylabel('z')
ax[1].set_ylabel('z')
ax[0].set_xlabel('Days')
ax[1].set_xlabel('Days')
#%%
plt.figure()
tl = range(0, 200)
plt.plot(um[tl,35:90:5], vm[tl,35:90:5])
plt.plot(ue[tl,35:90:5], ve[tl,35:90:5], linestyle='dashed')
plt.grid()
plt.xlabel('u')
plt.ylabel('v')
plt.axis('equal')
plt.tight_layout()
#%%
ztop = 30
upm = np.zeros((1420,z[ztop:-1].size))
for i in range(0, 1420):
 upm[i,:] = um[i,ztop:-1] - np.trapz(um[i,ztop:-1], x=z[ztop:-1], axis=-1)/(z[-1]-z[ztop])

deff = 1/2*time*np.trapz(upm[0:1407,:]**2, x=z[ztop:-1], axis=-1)/(z[-1]-z[ztop])
#%% ADVECTION DISTANCE
xe = integrate.cumtrapz(ue[0:1406,:], x=time[0:1406], axis=0)
ztop = 40#49
xem = integrate.trapz(xe[:,ztop:-1], x=z[ztop:-1], axis=-1)/(z[-1]-z[ztop])
xp2 = np.zeros((1405,))
for i in range(0, 1405):
    xp2[i] = integrate.trapz((xe[i,ztop:-1] - xem[i]), x=z[ztop:-1])**2
    
kappa2 = np.gradient(xp2)/np.gradient(time[0:1405])

plt.figure()
plt.plot(time/86400, xposunwrap[0:1407]-xposunwrap[0])
plt.plot(time[0:1405]/86400, xem)
plt.figure()
plt.plot(time/86400, XVar[0:1407] - XVar[0])
plt.plot(time[0:1405]/86400, xp2)
plt.plot(time[0:150]/86400, 0.5e1*(time[0:150]-time[0])**1)

plt.figure()
plt.plot(time/86400, kappa)
plt.plot(time[0:1405]/86400, kappa2)
plt.plot(time/86400, deff/5)

#plt.figure()
#plt.pcolor(time[0:1405]/86400,z[0:-2], np.transpose(xe[:,0:-2]), vmin=0, vmax=5e3)