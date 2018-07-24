#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 14:25:56 2018

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

filename = '/home/jacob/dedalus/LATMIX/front.mat'
f = h5py.File(filename, 'r')

# List all groups
print("Keys: %s" % f.keys())
keys = list(f.keys())[:]

vw = f['vw']
uw = f['uw']
wb = f['wflux']
u = f['u']
#v = f['v']
vm = f['va']
tx = f['tx']
ty = f['ty']
z = f['z'][:,0]
time = f['t'][:,0]
um = np.mean(u, axis=-1)
vwz = np.gradient(vw, axis=-1)/np.gradient(z)
uwz = np.gradient(uw, axis=-1)/np.gradient(z)
wbz = np.gradient(wb, axis=-1)/np.gradient(z)
q = f['q']
um = f['ua']
vm = f['va']
z = f['z'][:,0]
Bo = f['Bo']
vz = np.gradient(vm, axis=-1)/np.gradient(z)
uz = np.gradient(um, axis=-1)/np.gradient(z)
vt = np.transpose( np.gradient(np.transpose(vm), axis=-1)/np.gradient(time))
ut = np.transpose( np.gradient(np.transpose(um[0:time.size]), axis=-1)/np.gradient(time))
coriolis = 9.3e-5

#%%
tl = range(0,500)
plt.figure()
#plt.pcolor(time[0:300]/86400, z[0:-2], np.transpose(q[0:300,0:-2]), vmin=-1e-8, vmax=1e-8)
#plt.pcolor(time[tl]/86400, z[0:-2], np.transpose(vwz[tl,0:-2]), vmin=-1e-4/10, vmax=1e-4/10)
#plt.pcolor(time[0:300]/86400, z[0:-2], np.transpose(um[0:300,0:-2]), vmin=-0.2, vmax=0.2)
plt.pcolor(time[tl]/86400, z[0:-2], np.transpose(vt[tl,0:-2]), vmin=-1e-5, vmax=1e-5)

plt.colorbar()

#%% 
cm = 1e-5
tl = range(0,500)
summed = ut[:,:] - coriolis*vm[:,:] + uwz[:,:]
fig, ax = plt.subplots(4,2, sharex=True)
ax[0,0].pcolor(time[tl]/86400, z[0:-1], np.transpose(ut[tl,0:-1]), vmin=-cm, vmax = cm)
ax[1,0].pcolor(time[tl]/86400, z[0:-1], np.transpose(coriolis*vm[tl,0:-1]), vmin=-cm, vmax = cm)
ax[2,0].pcolor(time[tl]/86400, z[0:-1], np.transpose(-uwz[tl,0:-1]), vmin=-cm, vmax = cm)
ax[3,0].pcolor(time[tl]/86400, z[0:-1], np.transpose(summed[tl,0:-1]), vmin=-cm, vmax = cm)

# V 
summed = vt[:,:] + coriolis*um[0:time.size,:] + vwz[:,:]
ax[0,1].pcolor(time[tl]/86400, z[0:-1], np.transpose(vt[tl,0:-1]), vmin=-cm, vmax = cm)
ax[1,1].pcolor(time[tl]/86400, z[0:-1], np.transpose(-coriolis*um[tl,0:-1]), vmin=-cm, vmax = cm)
ax[2,1].pcolor(time[tl]/86400, z[0:-1], np.transpose(-vwz[tl,0:-1]), vmin=-cm, vmax = cm)
ax[3,1].pcolor(time[tl]/86400, z[0:-1], np.transpose(summed[tl,0:-1]), vmin=-cm, vmax = cm)

ax[3,0].set_xlabel('Days')
ax[3,1].set_xlabel('Days')

ax[0,0].set_ylabel('$u_t$')
ax[1,0].set_ylabel('$fv$')
ax[2,0].set_ylabel(r'$-\langle uw \rangle_z$')
ax[3,0].set_ylabel('Res')

ax[0,1].set_ylabel('$v_t$')
ax[1,1].set_ylabel('$-fu$')
ax[2,1].set_ylabel(r'$-\langle vw \rangle_z$')
ax[3,1].set_ylabel('Res')

plt.tight_layout()

#%%SHEAR BUDGET
vzt = np.transpose( np.gradient(np.transpose(vz), axis=-1)/np.gradient(time))
uzt = np.transpose( np.gradient(np.transpose(uz[0:time.size]), axis=-1)/np.gradient(time))
vwzz = np.gradient(vwz, axis=-1)/np.gradient(z)
uwzz = np.gradient(uwz, axis=-1)/np.gradient(z)
#%%
cm = 1e-6
cm = 5e-7
tl = range(0,500)
summed = uzt[:,:] - coriolis*vz[:,:] + uwzz[:,:]
fig, ax = plt.subplots(4,2, sharex=True)
ax[0,0].pcolor(time[tl]/86400, z[0:-2], np.transpose(uzt[tl,0:-2]), vmin=-cm, vmax = cm)
ax[1,0].pcolor(time[tl]/86400, z[0:-2], np.transpose(coriolis*vz[tl,0:-2]), vmin=-cm, vmax = cm)
ax[2,0].pcolor(time[tl]/86400, z[0:-2], np.transpose(-uwzz[tl,0:-2]), vmin=-cm, vmax = cm)
ax[3,0].pcolor(time[tl]/86400, z[0:-2], np.transpose(summed[tl,0:-2]), vmin=-cm, vmax = cm)

# V 
summed = vzt[:,:] + coriolis*uz[0:time.size,:] + vwzz[:,:]
ax[0,1].pcolor(time[tl]/86400, z[0:-2], np.transpose(vzt[tl,0:-2]), vmin=-cm, vmax = cm)
ax[1,1].pcolor(time[tl]/86400, z[0:-2], np.transpose(-coriolis*uz[tl,0:-2]), vmin=-cm, vmax = cm)
ax[2,1].pcolor(time[tl]/86400, z[0:-2], np.transpose(-vwzz[tl,0:-2]), vmin=-cm, vmax = cm)
ax[3,1].pcolor(time[tl]/86400, z[0:-2], np.transpose(summed[tl,0:-2]), vmin=-cm, vmax = cm)

ax[3,0].set_xlabel('Days')
ax[3,1].set_xlabel('Days')

ax[0,0].set_ylabel('$u_{zt}$')
ax[1,0].set_ylabel('$fv_z$')
ax[2,0].set_ylabel(r'$-\langle uw \rangle_{zz}$')
ax[3,0].set_ylabel('Res')

ax[0,1].set_ylabel('$v_{zt}$')
ax[1,1].set_ylabel('$-fu_z$')
ax[2,1].set_ylabel(r'$-\langle vw \rangle_{zz}$')
ax[3,1].set_ylabel('Res')

plt.tight_layout()

#%%
Bo = f['Bo']
EBF = integrate.trapz(um[:,70:98]*5e-7, x=z[70:98], axis=-1)
EBF = -ty[:]/coriolis*-5e-7
EBF = EBF[:,0]
tl = range(0, 400)
fig, ax = plt.subplots(2,1)
ax[0].plot(time[tl]/86400, tx[tl], label=r'$\tau_x$')
ax[0].plot(time[tl]/86400, ty[tl], label=r'$\tau_y$')
ax[0].legend()
ax[1].plot(time[tl]/86400, Bo[tl], label=r'$B_o$')
#ax[1].plot(time[tl]/86400, EBF[tl], label=r'$EBF$')

ax[1].legend()
ax[1].set_xlabel('Days')

#%% WIND WORK
tl = range(0,400)
WW = 0*tx[:,0]*um[0:tx.size,98] + ty[:,0]*vm[0:tx.size,98]
plt.figure()
plt.plot(time[tl]/86400, WW[tl])
WW = 0*tx[:,0]*umno[0:tx.size,98] + ty[:,0]*vmno[0:tx.size,98]

plt.plot(time[tl]/86400, WW[tl], linestyle='dashed')
plt.grid()
#%%
timen = time/(2*np.pi/coriolis)
Bo = f['Bo']
EBF = integrate.trapz(um[:,70:98]*5e-7, x=z[70:98], axis=-1)
EBF = -ty[:]/coriolis*-5e-7
EBF = EBF[:,0]
tl = range(0, 400)
fig, ax = plt.subplots(2,1)
ax[0].plot(timen[tl], tx[tl], label=r'$\tau_x$')

ax[0].plot(timen[tl], ty[tl], label=r'$\tau_y$')
ax[0].legend()
ax[1].plot(timen[tl], Bo[tl], label=r'$B_o$')
#ax[1].plot(time[tl]/86400, EBF[tl], label=r'$EBF$')

ax[1].legend()
ax[1].set_xlabel('Inertial Periods')
#%%
vtot = vm + -5e-7/coriolis*(z-z[1])
uv = uzno
vv = vzno
tl = range(0, 300)
plt.figure()
for i in range(50, 96, 10 ):
    plt.plot(uv[tl, i], vv[tl, i], label=str(z[i]), linewidth=2)
    j = np.where((time/86400>0.1))[0][0]
    plt.plot(uv[j,i], vv[j,i], marker='x')
#plt.plot(-0.1*np.cos(coriolis*time[tl])+0.0, 0.1*np.sin(coriolis*time[tl])+0.05 )
plt.axis('equal')
plt.grid()
plt.legend()

#%%
#VGZ = 5e-7
#vwzz = np.gradient(vwz, axis=-1)/np.gradient(z)
plt.figure()
plt.pcolor(time[tl]/86400, z[0:-2], np.transpose(-(vw[tl,0:-2]-0*uw[tl,0:-2])), vmin=-5e-4, vmax=5e-4)
plt.colorbar()
#plt.clim

#%%
zl = np.where((z>-50) & (z<-5))[0]
VGZ = 5e-7/9.3e-5
plt.figure()
#plt.plot(time/86400, integrate.trapz(vw[:,zl], x=z[zl], axis=-1)/(z[zl[-1]] -z[zl[0]])*VGZ, label=')
plt.plot(time/86400, np.max(vw[:,zl], axis=-1)*VGZ, label='Max(GSP)')
plt.plot(time/86400, -Bo[:,0]-EBF, label='$B_o$ + EBF')
plt.legend()

#%%
tl = range(0, 300)
cm = 1e-7
plt.figure()
plt.pcolor(time[tl]/86400,z, np.transpose(0*wbz[tl,:]-vwz[tl,:]*(-5e-7/coriolis)), vmin = -cm, vmax=cm)
plt.colorbar()
