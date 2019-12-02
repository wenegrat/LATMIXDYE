#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 14:33:35 2019

@author: jacob
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script makes the dye dispersal plots for figure 2 of the manuscript
(ie. dye dispersal and averaged vel profiles for the 3 different runs)

Created on Wed Mar 13 10:05:25 2019

@author: jacob
"""
import h5py
import scipy.io as spio
import matplotlib.pyplot as plt
import numpy as np
from cmocean import cm as cmo
import scipy.integrate as integrate 
import matplotlib.ticker 
import scipy.interpolate as interpolate

plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 20
plt.rcParams['contour.negative_linestyle'] = 'solid'
#matplotlib.rcParams.update({'font.size': 22})
#%% LOAD DATA

filename = '/home/jacob/dedalus/LATMIX/run10k.mat' # Front run
f = h5py.File(filename, 'r')

# List all groups
print("Keys: %s" % f.keys())
keys = list(f.keys())[:]

u_f = f['u']
va_f = f['va']

dye1_f = f['dye1']
time_f = f['t'][:,0]
b_f = f['b'] # total buoyancy
x = f['x'][:,0]
z = f['z'][:,0]
nt, nz, nx = dye1_f.shape


#%% LOAD MET DATA

filename = '/home/jacob/dedalus/LATMIX/LES_COARE35_forcing.mat' # Front run
matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)

# List all groups
yearday = matfile['yearday']
tau_cs = matfile['tau_cs']
tau_d = matfile['tau_d']
qnet = matfile['net_heat_flux']
tmag = np.sqrt(tau_cs**2 + tau_d**2)


#%% INTERPOLATE ALONG TIME DIMENSION (buoyancy and dye1)
timesi = np.linspace(0, 66-64.87, 1000)
nti = timesi.size
dyei = np.zeros((nti, nz, nx))
bi = np.zeros((nti, nz, nx))
for i in range(0, nx):
    print(i/nx)

    f = interpolate.interp1d(time_f/86400+64.5-64.87, dye1_f[0:686,:,i], axis=0)
    dyei[:,:,i] = f(timesi)
    f = interpolate.interp1d(time_f/86400+64.5-64.87, b_f[0:686,:,i], axis=0)
    bi[:,:,i] = f(timesi)
#%%


cmap = 'gnuplot'
save = True
conts = np.linspace(-3, 0, 50)

bconts = np.linspace(-0.005, 0.02, 25)
btemp = b_f[162,74,1280] # buoyancy at center of initial dye patch
bind = np.argmin(np.abs(bconts-btemp))
bconts = bconts + (bconts[bind] - btemp)# + 0.5*(bconts[1]-bconts[0])
fac = 1.5
plt.rcParams['font.size'] = 20/fac

#indb = np.argmin(np.abs(time_f/86400-64.87+64.5))
for i in range(0, nti, 1):
    print(i)
    fig = plt.figure(figsize=(16/fac,9/fac), dpi=fac*100)
    axdye = plt.subplot2grid((3,3), (0,0), rowspan = 2, colspan=3)

#    dyevar = np.log10(dye1_f[i,1:,:])
    dyevar = np.log10(dyei[i,1:,:])
    dyevar[np.isnan(dyevar)] = -10
    dyevar[~np.isfinite(dyevar)] = -10
    
#    bvar = b_f[i,:-1,:]
    bvar = bi[i,:-1,:]
    im = axdye.contourf(x/1e3, z[:-1], dyevar, conts, extend='min', cmap = cmap)
    for c in im.collections:
        c.set_edgecolor("face")
    axdye.contour(x/1e3, z[:-1], bvar, bconts, colors='w',  linewidths=0.5)
    axdye.set_xlabel('y [km]')
    axdye.set_ylabel('z [m]')
    axdye.set_yticks([-150, -100, -50, 0])
    axdye.set_xticks([0, 2, 4, 6, 8, 10])
    axdye.set_ylim(-150, 0)
#    axdye.set_title(f'Days since dye release: {time_f[i]/86400 - time_f[indb]/86400:3.2f}')
    axdye.set_title(f'Days since dye release: {timesi[i]:3.2f}')

    axflux = plt.subplot2grid((3,3), (2,0), colspan=3)
    color = 'tab:blue'
    axflux.plot(yearday, tmag, color=color, linewidth=2)
    axflux.tick_params(axis='y', labelcolor=color)
    axflux.set_ylabel('Wind-stress\n magnitude\n [N m$^{-2}$]', color=color)  # we already handled the x-label with ax1
    axflux.set_ylim((0, 1.2))
    axflux.set_yticks([0, 0.3, 0.6, 0.9,  1.2])
    axflux.set_xlim((64.5, 66))
    axflux.set_xlabel('yearday')
    
    axflux2 = axflux.twinx()
#    axflux.axvline(time_f[i]/86400+64.5, linestyle='dashed', color='k')
    axflux.axvline(timesi[i]+64.87, linestyle='dashed', color='k')

    color = 'tab:red'
    axflux2.plot(yearday, qnet, color=color, linewidth=2)
    axflux2.set_ylabel('Net surface\n heat flux\n [W m$^{-2}$]', color=color)  # we already handled the x-label with ax1
    axflux2.tick_params(axis='y', labelcolor=color)
    axflux2.set_ylim((0, 1200))
    axflux2.set_xlim((64.87, 66))
    axflux2.set_yticks([0, 300, 600, 900,  1200])
    axflux.grid(linestyle='dashed')
#    axflux.text( -0.22, 0.75 ,'c)', color='k', size=20, bbox=dict(facecolor='w', edgecolor='k'),transform=axflux.transAxes )

    plt.subplots_adjust(wspace=0, hspace=0.5)
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.91, 0.45, 0.01, 0.45])
    cb = fig.colorbar(im, cax=cbar_ax, extend='min')
    cb.set_ticks([-3, -2, -1, 0])
    cb.set_label('log$_{10}$[Normalized Concentration]')
    fig.tight_layout(pad=0.2)
    if save:
        plt.savefig(f'/home/jacob/Dropbox/GulfStreamDye/MovieFigs/LatmixDyeModel_{i}.png')
        plt.close('all')
    
## TO MAKE MOVIE
# ffmpeg -start_number 0 -i LatmixDyeModel_%d.png -pix_fmt yuv420p -crf 30 dye_model.mov
