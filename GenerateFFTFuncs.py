#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 15:32:30 2018

@author: jacob
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.interpolate as interpolate
#%%
moviefile = np.loadtxt('/home/jacob/LATMIX/vapor/dyeflux.txt')

NX=1280
NY=100
N_TH=0
LX=5000
LY=150



# Get the number of timesteps (note this depends on how the array is shaped...)
nk=(moviefile.shape[0])/(NX*NY);
nk = np.int(nk)
count=0;

# zero matrix mat
uth = np.zeros((NX, NY, nk))
dthdx = np.zeros((NX,NY, nk))
wth = np.zeros((NX, NY, nk))
dthdz = np.zeros((NX,NY, nk))
dtime = np.zeros((NX,NY, nk))

# Make sure this order matches how it is being saved in the .f file
for k in range(0, nk):
    print(k)
    for j in range(0, NY):
        for i in range(0, NX):
            dtime[i,j,k]=moviefile[count,0]
            uth[i,j,k]=moviefile[count,1]
            dthdx[i,j,k] = moviefile[count,2] - 5e-7
            wth[i,j,k]=moviefile[count,3]
            dthdz[i,j,k] = moviefile[count,4]
            count=count+1

dtime = dtime[0,0,:]
#%%
kappah = -uth/dthdx
kappah[np.isinf(kappah)] = 0

#kappahx = -np.mean(uth, axis=1)/np.mean(dthdx, axis=1)
#kappahz = -np.mean(uth, axis=0)/np.mean(dthdx, axis=0)
#
#kappahx[np.isinf(kappahx)] = 0
#kappaz = -np.mean(wth[:,:,20:40],axis=-1)/np.mean(dthdz[:,:,20:40],axis=-1)
#kappaz = -wth/dthdz
zl = range(60, 98)
kappaz = -(np.mean(1/np.abs(z[zl[0]])*integrate.trapz(wth[:,zl,:], x=z[zl], axis=1), axis=0)
            /np.mean(1/np.abs(z[zl[0]])*integrate.trapz(dthdz[:,zl,:], x=z[zl], axis=1), axis=0))

kappah = -(np.mean(1/np.abs(z[zl[0]])*integrate.trapz(uth[:,zl,:], x=z[zl], axis=1), axis=0)
            /np.mean(1/np.abs(z[zl[0]])*integrate.trapz(dthdx[:,zl,:], x=z[zl], axis=1), axis=0))

#ts = range(40, 50)
#kappaz = - np.mean(np.mean(wth[:,:,ts], axis=-1), axis=0)/np.mean(np.mean(dthdz[:,:,ts], axis=-1), axis=0)
kappaz[np.isinf(kappaz)] = 0

#kappa
#%%
def running_mean(x, N):
#    cumsum = np.cumsum(np.insert(x, 0, 0), axis=-1) 
    cumsum = np.cumsum(x, axis=-1) 

    return (cumsum[:,N:] - cumsum[:,:-N]) / float(N)
zl = range(25, 98)
tsm = 15
smf = running_mean(1/np.abs(z[zl[0]])*integrate.trapz(uth[:,zl,:], x=z[zl], axis=1), tsm)
smg = running_mean(1/np.abs(z[zl[0]])*integrate.trapz(dthdx[:,zl,:], x=z[zl], axis=1), tsm)

smf[-1,:] = 0
smg[-1,:] = 0
kaptest = -smf/(smg)
kaptest[np.isinf(kaptest)] = 0
kaptest[np.isnan(kaptest)] = 0
#%%
plt.figure()
plt.pcolor( np.transpose(uth[:,:,2]))
plt.colorbar()
plt.clim((-1e-3, 1e-3))

#%%
#zl = range(60, 98)

plt.figure()
plt.plot(time, kappa)
plt.plot(dtime, kappah, marker='x')