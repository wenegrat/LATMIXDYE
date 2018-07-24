#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:50:08 2018

@author: jacob
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
import scipy.integrate as integrate
import scipy.interpolate as interpolate
#%%
moviefile = np.loadtxt('/home/jacob/LATMIX/vapor/eterms.txt')

NX=1280/3
NX = np.ceil(NX)
NX = NX.astype('int')
NY=100
N_TH=0
LX=5000
LY=150



# Get the number of timesteps (note this depends on how the array is shaped...)
nk=(moviefile.shape[0])/(NX*NY);
nk = np.int(nk)
count=0;

# zero matrix mat
vw = np.zeros((NX, NY, nk))
uu = np.zeros((NX, NY, nk))
uw = np.zeros((NX, NY, nk))
ww = np.zeros((NX, NY, nk))
wd = np.zeros((NX, NY, nk))
kx = np.zeros((NX, NY, nk))
wb = np.zeros((NX, NY, nk))
dtime = np.zeros((NX,NY, nk))

# Make sure this order matches how it is being saved in the .f file
for k in range(0, nk):
    print(k)
    for j in range(0, NY):
        for i in range(0, NX):
            dtime[i,j,k]=moviefile[count,0]
            kx[i,j,k]=moviefile[count,1]
            vw[i,j,k]=moviefile[count,2]
            wb[i,j,k]=moviefile[count,3]
            uu[i,j,k]=moviefile[count,4]
            uw[i,j,k]=moviefile[count,5]
            ww[i,j,k]=moviefile[count,6]
            wd[i,j,k]=moviefile[count,7]
            count=count+1

k = kx[:,0,0]
dtime = dtime[0,0,:]
#%%
filename = '/home/jacob/dedalus/LATMIX/nofront.mat'
f = h5py.File(filename, 'r')
vm = f['va']
z  = f['z'][:,0]
time = f['t'][:,0]
vz = np.gradient(vm, axis=-1)/np.gradient(z)
vza = vz - 0*(-5e-7/9.3e-5)
vzai = interpolate.interp1d(time, vza, axis=0)
vza = vzai(dtime[0:-2])
#%%

plt.rcParams.update({'font.size': 22})
plt.figure()
zl = range(40,78)
zl = range(50, 93)
tl = range(0, 14)
tl = range(30, 35)
#tl = range( 0, 4)
vwt = integrate.trapz(vw[:,:,tl], x=dtime[tl], axis=-1)/(dtime[tl[-1]]-dtime[tl[0]])
wbt = integrate.trapz(wb[:,:,tl], x=dtime[tl], axis=-1)/(dtime[tl[-1]]-dtime[tl[0]])
vwi = integrate.trapz(vwt[:,zl]*(5e-7/9.3e-5), x=z[zl], axis=1)/(z[zl[-1]] - z[zl[0]])
wbi = integrate.trapz(wbt[:,zl], x=z[zl], axis=1)/(z[zl[-1]] - z[zl[0]])
#plt.semilogx(k,uwi)
plt.semilogx(k[0:-1], integrate.cumtrapz(vwi[::-1])[::-1],
             label='$-S^2/f\int_{k_z}\int_{k_x}^{k_x^{Max}}\hat{v}^*\hat{w}\mathrm{dk_xdk_z}$') # Note integration begins at high wavenumbers
#plt.semilogx(k[0:], 2*np.pi/5e3*np.cumsum(uwi[::-1], axis=0)[::-1],
#             label='$-S^2/f\int_{k_z}\int_{k_x^{Max}}^{k_x}\hat{u}^*\hat{w}\mathrm{dk_xdk_z}$') # Note integration begins at high wavenumbers

plt.semilogx(k[0:-1], integrate.cumtrapz(wbi[::-1])[::-1], 
              label='$\int_{k_z}\int_{k_x}^{k_x^{Max}}\hat{w}^*\hat{b}\mathrm{dk_xdk_z}$') # Note integration begins at high wavenumbers
plt.xlabel('$k_x$')
plt.xlim((1e-3, k[-1]))
#plt.ylim((-3e-15, 3e-15))
plt.ylabel('$W/kg$')
plt.legend()
plt.grid()
plt.tight_layout()

#%% HOVMOLLER
zl = range(80, 95)
zl = range(40, 80)
vwi = integrate.trapz(vw[:,zl,:]*(5e-7/9.3e-5), x=z[zl], axis=1)/(z[zl[-1]] - z[zl[0]])
vwi = integrate.cumtrapz(vwi[::-1,:], axis=0)[::-1,:]
wbi = integrate.trapz(wb[:,zl,:], x=z[zl], axis=1)/(z[zl[-1]] - z[zl[0]])
wbi = integrate.cumtrapz(wbi[::-1,:], axis=0)[::-1,:]

vwin = np.zeros((k.size-1, dtime.size))
wbin = np.zeros((k.size-1, dtime.size))
for i in range(0,dtime.size):
    vwin[:,i] = vwi[:,i]/np.max(np.abs(vwi[:,i])) #Normalize by maximum integrated GSP
    wbin[:,i] = wbi[:,i]/np.max(np.abs(wbi[:,i]))

fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6))

ax[0].pcolor(dtime/86400, k, vwin, vmin=-1, vmax=1)
ax[0].set_yscale('log')
ax[0].set_ylim((1e-3, k[-1]))
ax[0].set_xlim((dtime[0]/86400, dtime[-1]/86400))
ax[1].pcolor(dtime/86400, k, wbin, vmin=-1, vmax=1)
ax[1].set_yscale('log')
ax[1].set_ylim((1e-3, k[-1]))
ax[1].set_xlim((dtime[0]/86400, dtime[-1]/86400))
ax[1].set_xlabel('Days')
ax[1].set_ylabel('$k_x$ (rad/m)')
ax[0].set_ylabel('$k_x$ (rad/m)')

ax[1].set_title('Normalized $\int_{k_x}^{k_x^{Max}}BFLUX\mathrm{dk_x}$')
ax[0].set_title('Normalized $\int_{k_x}^{k_x^{Max}}GSP\mathrm{dk_x}$')
plt.tight_layout()
#%% HOVMOLLER NOT CUMULATIVE
zl = range(80, 95)
zl = np.where((z>-75) & (z<-5))[0]
#zl = range(40, 80)
vwi = integrate.trapz(vw[:,zl,:]*(5e-7/9.3e-5), x=z[zl], axis=1)/(z[zl[-1]] - z[zl[0]])
#vwi = integrate.trapz(vw[:,zl,0:-2]*np.transpose(vza[:,zl]), x=z[zl], axis=1)/(z[zl[-1]] - z[zl[0]])

wbi = integrate.trapz(wb[:,zl,:], x=z[zl], axis=1)/(z[zl[-1]] - z[zl[0]])

vwin = np.zeros((k.size, dtime.size))
wbin = np.zeros((k.size, dtime.size))
for i in range(0,dtime.size-2):
    vwin[:,i] = k*vwi[:,i]#/np.max(np.abs(k*vwi[:,i])) #Normalize by maximum integrated GSP
    wbin[:,i] = k*wbi[:,i]#/np.max(np.abs(k*wbi[:,i]))

fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12,6))
tl = range(0, 20)
cl = 5e-10
ax[0].pcolor(dtime[tl]/86400, k, vwin[:,tl], vmin=-cl, vmax=cl)
ax[0].set_yscale('log')
ax[0].set_ylim((1e-3, k[-1]))
ax[0].set_xlim((dtime[0]/86400, dtime[-1]/86400))
ax[1].pcolor(dtime[tl]/86400, k, wbin[:,tl])
ax[1].set_yscale('log')
ax[1].set_ylim((1e-3, k[-1]))   
ax[1].set_xlim((dtime[tl[0]]/86400, dtime[tl[-1]]/86400))
ax[1].set_xlabel('Days')
ax[1].set_ylabel('$k_x$ (rad/m)')
ax[0].set_ylabel('$k_x$ (rad/m)')

ax[1].set_title('$k_xBFLUX$')
ax[0].set_title('$k_xGSP$')
plt.tight_layout()

#%% HOVMOLLER UU
zl = range(80, 95)
zl = np.where((z>-75) & (z<-5))[0]
vari = uu/np.max(np.abs(uu),axis=0)# 
#vari = wd/np.max(np.abs(wb), axis=0)
#vari= vw/np.max(np.abs(vw), axis=0)#
vari = uw**2/np.sqrt(uu*ww)
uui = integrate.trapz(vari[:,zl,:], x=z[zl], axis=1)/(z[zl[-1]] - z[zl[0]])

uuin = np.zeros((k.size, dtime.size))
for i in range(0,dtime.size):
    uuin[:,i] = uui[:,i]/np.max(np.abs(uui[:,i])) #Normalize by maximum integrated GSP

plt.figure(figsize=(12,3))

plt.pcolor(dtime/86400, k, uuin, vmin=-1, vmax=1)
plt.yscale('log')
plt.ylim((1e-3, k[-1]))
plt.xlim((dtime[0]/86400, dtime[-1]/86400))

plt.xlim((dtime[0]/86400, dtime[-1]/86400))
plt.xlabel('Days')
plt.ylabel('$k_x$ (rad/m)')
plt.colorbar()
#plt.title('$k_x(uu)$')
plt.tight_layout()
#%% DEPTH WAVENUMBER SLICE
tl = range(0, 14)
tl = range(30, 40)
tl = np.where((dtime/86400>1.) & (dtime/86400<1.5) )[0]
#tl = range( 0, 4)
vwt = integrate.trapz(vw[:,:,tl], x=dtime[tl], axis=-1)/(dtime[tl[-1]]-dtime[tl[0]])
wbt = integrate.trapz(wb[:,:,tl], x=dtime[tl], axis=-1)/(dtime[tl[-1]]-dtime[tl[0]])

for i in range(0, z.size):
    vwt[:,i] = vwt[:,i]*k
    wbt[:,i] = wbt[:,i]*k
    
fig, ax = plt.subplots(2,1, sharex=True)

ax[0].pcolor(k, z, np.transpose(vwt))
ax[0].set_xscale('log')
ax[0].set_xlim((1e-3, k[-1]))
ax[0].set_title('GSP, Days: %1.1f' %(dtime[tl[1]]/86400)  + ' - %1.1f'%(dtime[tl[-1]]/86400) )
ax[0].plot(k, -2*np.pi/k)
ax[0].set_ylim((z[0], 0))
ax[1].pcolor(k, z, np.transpose(wbt))
ax[1].set_xscale('log')
ax[1].set_xlim((1e-3, k[-1]))
ax[1].set_xlabel('$k_x$ (rad/m)')
ax[1].set_ylabel('z (m)')
ax[0].set_ylabel('z (m)')
ax[1].set_title('BFLUX')
plt.tight_layout()
#%%
#plt.figure()
#zl = range(40, 78)
#tl = range(0, 14)
#uwi = np.mean(integrate.trapz(-uw[:,zl,0:nk-1]*(5e-7/1e-4),x=z[zl], axis=1), axis=-1)
#wbi = np.mean(integrate.trapz(wb[:,zl,0:nk-1], x=z[zl], axis=1), axis=-1)
##plt.semilogx(k,uwi)
#plt.semilogx(k[0:], uwi) # Note integration begins at high wavenumbers
#plt.semilogx(k[0:], wbi) # Note integration begins at high wavenumbers

