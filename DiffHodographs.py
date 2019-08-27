#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 15:16:58 2018
Make hodographs with and without thermal wind mixing.

@author: jacob
"""
import h5py
import scipy.io as spio
import matplotlib.pyplot as plt
import numpy as np
from cmocean import cm as cmo

plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 14
plt.rcParams['contour.negative_linestyle'] = 'solid'
#%% LOAD NL SIM SNAPS
filename = '/home/jacob/dedalus/LATMIX/notw.mat'
fno = h5py.File(filename, 'r')
z = fno['z'][:,0]
umno = fno['ua']
vmno = fno['va']
vwno = fno['vw']
uwno = fno['uw']
qno = fno['q']
uzno = np.gradient(umno, axis=-1)/np.gradient(z)
vzno = np.gradient(vmno, axis=-1)/np.gradient(z)
timeno=fno['t'][:,0]
#qzno = fno['qz']
filename = '/home/jacob/dedalus/LATMIX/run10k.mat'
f = h5py.File(filename, 'r')

um = f['ua']
vm = f['va']
vw = f['vw']
uw = f['uw']
q = f['q']
#qz = f['qz']
#ud = um[:,:] - umno[0:1569,:]
#vd = vm[:,:] - vmno[0:1569,:] 
#
#udz = np.gradient(ud, axis=-1)/np.gradient(z)
#vdz = np.gradient(vd, axis=-1)/np.gradient(z)

vwz = np.gradient(vw, axis=-1)/np.gradient(z)
uwz = np.gradient(uw, axis=-1)/np.gradient(z)
vwzno = np.gradient(vwno, axis=-1)/np.gradient(z)
uwzno = np.gradient(uwno, axis=-1)/np.gradient(z)


uz = np.gradient(um, axis=-1)/np.gradient(z)
vz = np.gradient(vm, axis=-1)/np.gradient(z)
time=f['t'][:,0]
z = f['z'][:,0]

vt = np.transpose( np.gradient(np.transpose(vm[0:time.size]), axis=-1)/np.gradient(time))
ut = np.transpose( np.gradient(np.transpose(um[0:time.size]), axis=-1)/np.gradient(time))
vtno = np.transpose( np.gradient(np.transpose(vmno[0:time.size]), axis=-1)/np.gradient(time))
utno = np.transpose( np.gradient(np.transpose(umno[0:time.size]), axis=-1)/np.gradient(time))

utd = ut - utno
vtd = vt - vtno
#
#uwzd = uwz - uwzno[0:time.size,:]
#vwzd = vwz - vwzno[0:time.size,:]

coriolis = 9.3e-5
dbdx = -5e-7
Vgz = dbdx/coriolis
Vg = (z[1:]-z[1])*Vgz

#%% SINGLE DEPTH

tl = range(0,500)
timen = time[tl]/(2*np.pi/coriolis)
zdepth = 85
fig,ax = plt.subplots(2,1, sharex=True)
ax[0].plot(timen, ut[tl,zdepth], label='d/dt')
ax[0].plot(timen, coriolis*vm[tl,zdepth], label='Cori')
ax[0].plot(timen, -uwz[tl,zdepth], label='Stress')
ax[0].legend()

ax[1].plot(timen, vt[tl,zdepth], label='d/dt')
ax[1].plot(timen, -coriolis*um[tl,zdepth], label='Cori')
ax[1].plot(timen, -vwz[tl,zdepth], label='Stress')
ax[1].legend()

ax[1].set_xlim((0, 1))
ax[0].set_xlim((0, 1))
ax[0].set_title(z[zdepth])
#%% HODOGRAPHS FOR SUPPLEMENTARY INFORMATION
yls = (-0.2, 0.2)
fog, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(9, 6.5))
tl = range(0, 269)
tl = 0
tr = np.argmin(np.abs(time/86400 + 64.5 - 65.2))
trno = np.argmin(np.abs(timeno/86400 + 64.5 - 65.2))
for i in range(0, 7):
    zind = np.argmin(np.abs(z + i*10 + 10))
    ax1.plot(-vm[tl:tr, zind], um[tl:tr,zind])
    ax2.plot(-vmno[tl:trno, zind], umno[tl:trno,zind], label='{} m'.format(str(np.round(z[zind]/10)*10)))

#ax2.plot( -vmno[tl, 37:95:10]+0*Vg[37:95:10], umno[tl,37:95:10])
ax2.legend(loc='center left',bbox_to_anchor=(1, 0.5))
ax1.set_ylim(yls)
ax1.set_xlim(yls)
ax1.grid()
ax2.set_ylim(yls)
ax2.set_xlim(yls)
ax2.grid()
ax1.set_ylabel('Across front velocity [m/s]')
ax1.set_xlabel('Along front velocity [m/s]')
ax2.set_xlabel('Along front velocity [m/s]')
ax1.set_title('FRONT')
ax2.set_title('NO-GEOMIX')
ax1.set(adjustable='box-forced', aspect='equal')
ax2.set(adjustable='box-forced', aspect='equal')
plt.tight_layout()
#plt.axes('equal')
#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/FigureSI1.pdf', bbox_inches='tight')

#%%

tl = range(0, 300)
plt.figure()
for i in range(37, 95, 10 ):
    plt.plot(umno[tl, i], vmno[tl, i], label=str(z[i]), linewidth=2, linestyle='dashed')
    plt.plot(um[tl, i], vm[tl, i], label=str(z[i]), linewidth=2)

#plt.plot(-0.1*np.cos(coriolis*time[tl])+0.0, 0.1*np.sin(coriolis*time[tl])+0.05 )
plt.axis('equal')
plt.grid()
plt.legend()


#%%
tl = range(0, 210)
plt.figure()
for i in range(80, 98, 4 ):
#        plt.plot(uz[tl, i] ,vz[tl,i]-5e-7/9.3e-5, label=str(z[i]), linewidth=2)

    plt.plot(ud[tl, i] ,vd[tl,i], label=str(z[i]), linewidth=2)
    zd = 10
    j = np.where((time/86400>0.31))[0][0]
    plt.plot(ud[j,i], vd[j,i], marker='x')
#    plt.plot((ud[tl, i]-ud[tl,i-zd])/(z[i] - z[i-zd]), (vd[tl, i] - vd[tl, i-zd])/(z[i] - z[i-zd]) - 5e-7/9.3e-5, label=str(z[i]), linewidth=2)

#plt.plot(-0.1*np.cos(coriolis*time[tl])+0.0, 0.1*np.sin(coriolis*time[tl])+0.05 )
plt.axis('equal')
plt.grid()
plt.legend()


#%% MOM DIFF PLOT
#%% 


cm = 1e-5
tl = range(0,500)
summed = utd[:,:] - coriolis*vd[:,:] + uwzd[:,:]
timen = time[tl]/(2*np.pi/coriolis)
fig, ax = plt.subplots(4,2, sharex=True)
ax[0,0].pcolor(timen, z[0:-2], np.transpose(utd[tl,0:-2]), vmin=-cm, vmax = cm)
ax[1,0].pcolor(timen, z[0:-2], np.transpose(coriolis*vd[tl,0:-2]), vmin=-cm, vmax = cm)
ax[2,0].pcolor(timen, z[0:-2], np.transpose(-uwzd[tl,0:-2]), vmin=-cm, vmax = cm)
ax[3,0].pcolor(timen, z[0:-2], np.transpose(summed[tl,0:-2]), vmin=-cm, vmax = cm)

# V 
summed = vtd[:,:] + coriolis*ud[0:time.size,:] + vwzd[:,:]
ax[0,1].pcolor(timen, z[0:-2], np.transpose(vtd[tl,0:-2]), vmin=-cm, vmax = cm)
ax[1,1].pcolor(timen, z[0:-2], np.transpose(-coriolis*ud[tl,0:-2]), vmin=-cm, vmax = cm)
ax[2,1].pcolor(timen, z[0:-2], np.transpose(-vwzd[tl,0:-2]), vmin=-cm, vmax = cm)
ax[3,1].pcolor(timen, z[0:-2], np.transpose(summed[tl,0:-2]), vmin=-cm, vmax = cm)

ax[3,0].set_xlabel('Inertial Periods')
ax[3,1].set_xlabel('Inertial Periods')

ax[0,0].set_ylabel('$u_{t}$')
ax[1,0].set_ylabel('$fv$')
ax[2,0].set_ylabel(r'$-\langle uw \rangle_{z}$')
ax[3,0].set_ylabel('Res')

ax[0,1].set_ylabel('$v_{t}$')
ax[1,1].set_ylabel('$-fu$')
ax[2,1].set_ylabel(r'$-\langle vw \rangle_{z}$')
ax[3,1].set_ylabel('Res')

ax[0,0].set_xlim((0, 1))
plt.tight_layout()

#%%
plt.figure()
plt.plot(timen[tl], tm[tl])
#%% SINGLE DEPTH

tl = range(0,500)
timen = time[tl]/(2*np.pi/coriolis)
zdepth = 85
fig,ax = plt.subplots(2,1, sharex=True)
ax[0].plot(timen, utd[tl,zdepth], label='d/dt')
ax[0].plot(timen, coriolis*vd[tl,zdepth], label='Cori')
ax[0].plot(timen, -uwzd[tl,zdepth], label='Stress')
ax[0].legend()

ax[1].plot(timen, vtd[tl,zdepth], label='d/dt')
ax[1].plot(timen, -coriolis*ud[tl,zdepth], label='Cori')
ax[1].plot(timen, -vwzd[tl,zdepth], label='Stress')
ax[1].legend()

ax[1].set_xlim((0, 1))
ax[0].set_xlim((0, 1))
ax[0].set_title(z[zdepth])

#%%
cl = 5e-9
tl = range(0,500)
timen = time[tl]/(2*np.pi/coriolis)
fig, ax = plt.subplots(3,1, sharex = True)

ax[0].pcolor(timen, z[2:-2], np.transpose(q[0, tl, 2:-2]), vmin=-cl, vmax=cl)
ax[1].pcolor(timen, z[2:-2], np.transpose(qno[0, tl, 2:-2]), vmin=-cl, vmax=cl)
ax[2].pcolor(timen, z[2:-2], np.transpose(q[0, tl, 2:-2] - qno[0, tl,2:-2]), vmin=-cl, vmax=cl)

#%%
tl = range(0,500)
timen = time[tl]/(2*np.pi/coriolis)
zl = range(90, 98)
qi = integrate.trapz(q[0, :, zl], x=z[zl], axis=-1)/(z[zl[-1]] - z[zl[0]])
qnoi = integrate.trapz(qno[0,:,zl], x=z[zl], axis=-1)/(z[zl[-1]] - z[zl[0]])

qit = np.gradient(qi)/np.gradient(time)
qnoit = np.gradient(qnoi[0:time.size])/np.gradient(time)
plt.figure()
plt.plot(timen, qi[tl], label='Front')
plt.plot(timen, qnoi[tl], label='NoTW')
plt.xlabel('Inertial Periods')
plt.ylabel('PV (s$^{-3}$)')
plt.legend()
plt.tight_layout()