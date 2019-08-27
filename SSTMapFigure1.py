#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 15:51:49 2019

@author: jacob
"""

import h5py
import scipy.io as spio
#import mpl_toolkits.basemap 
import matplotlib.pyplot as plt
import numpy as np
from cmocean import cm as cmo
from scipy.interpolate import griddata
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import matplotlib.patches as patches

#from matplotlib.mlab import griddata
import datetime as dt
import gsw
import scipy.interpolate as interp
from matplotlib.gridspec import GridSpec

from ReadAtlantis_ALL import ReadAtlantis_ALL
from ReadAtlantis import ReadAtlantisSection
from ReadKnorr_ALL import ReadKnorr_ALL

#%% LOAD SST MAP

filename = '/home/jacob/dedalus/LATMIX/SatData/AVHRR493data.mat'
matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)

sstdatafile=matfile['sstdata']

temp0 = sstdatafile[0]
temp1 = sstdatafile[1]
temp2 = sstdatafile[2]

sstdata = temp0.sst
sst=(sstdata.bytes-0.5)*sstdata.scale; 
sst[sst<0] = np.nan; 
sst=sst+sstdata.min

latgrid = temp0.subset.latgrid
longrid = temp0.subset.lonEgrid

latlimr = np.argmin(np.abs(latgrid-38))
latliml = np.argmin(np.abs(latgrid-40))
lonliml = np.argmin(np.abs(longrid+66))
lonlimr = np.argmin(np.abs(longrid+62))

#%% LOAD MET DATA

filename = '/home/jacob/dedalus/LATMIX/LES_COARE35_forcing.mat' # Front run
matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)

# List all groups
yearday = matfile['yearday']
tau_cs = matfile['tau_cs']
tau_d = matfile['tau_d']
qnet = matfile['net_heat_flux']
tmag = np.sqrt(tau_cs**2 + tau_d**2)

#%% LOAD FLOAT DATA
filename = '/home/jacob/dedalus/LATMIX/FloatData/Mar05_SI_2_Track.mat'
matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)
floatstruct = matfile['F']
flat = floatstruct.lat
flon = floatstruct.lon
fyd = floatstruct.yd-1 #Craig convention

fyds = floatstruct.yds - 1
flats = floatstruct.lats
flons = floatstruct.lons
#%% LOAD KNARR TRACK


filename = '/home/jacob/dedalus/LATMIX/LatMix_2012/my_triaxus_SI2'
matfile = spio.loadmat(filename,struct_as_record=True, squeeze_me=True)

jday_ts = matfile['jday_ts']

latK = matfile['lat_ts']
lonK = matfile['lon_ts']

#%% LOAD ATLANTIS TRACK
latA, lonA, timeA = ReadAtlantis_ALL()

#%%

cols = 7

fig = plt.figure(figsize=(12.22, 6.16))

# SST PLOT

aspratio = 1/np.cos(39*np.pi/180)
axsst = plt.subplot2grid((3,cols), (0,0), rowspan = 2, colspan=4)
isst = axsst.pcolor(longrid[lonliml:lonlimr], latgrid[latliml:latlimr], np.transpose(sst[lonliml:lonlimr, latliml:latlimr]), cmap=cmap, vmin = tliml, vmax=tlimh)
axsst.plot(flons, flats, linewidth=4)
axsst.plot(lonK, latK, linewidth=4)
axsst.plot(lonA, latA, linewidth=4)
indrelease = np.argmin(np.abs(fyds - 64.86))
axsst.plot(flons[indrelease], flats[indrelease], marker='x')
#plt.colorbar(isst, cax = axsst)
axsst.set_xlim(-66, -62)
axsst.set_ylim(38, 40)
axsst.set_aspect(aspratio)
axsst.set_ylabel('Latitude')
axsst.set_xlabel('Longitude')
axsst.grid(linestyle='dashed')
# SURFACE FORCING PLOT
axflux = plt.subplot2grid((3,cols), (2,1), colspan=cols-2)
color = 'tab:blue'
axflux.plot(yearday, tmag, color=color, linewidth=2)
axflux.tick_params(axis='y', labelcolor=color)
axflux.set_ylabel('Wind-stress\n magnitude\n (N m$^{-2}$)', color=color)  # we already handled the x-label with ax1
axflux.set_ylim((0, 1.2))
axflux.set_yticks([0, 0.3, 0.6, 0.9,  1.2])
axflux.set_xlim((64.5, 66))
axflux.set_xlabel('yearday')

axflux2 = axflux.twinx()
axflux.axvline(64.86, linestyle='dashed', color='k')
color = 'tab:red'
axflux2.plot(yearday, qnet, color=color, linewidth=2)
axflux2.set_ylabel('Net surface\n heat flux\n (W m$^{-2}$)', color=color)  # we already handled the x-label with ax1
axflux2.tick_params(axis='y', labelcolor=color)
axflux2.set_ylim((0, 1200))
axflux2.set_yticks([0, 300, 600, 900,  1200])
axflux.grid(linestyle='dashed')

# SALINITY SECTION
axSal = plt.subplot2grid((3,cols), (0,4), colspan=cols-4, rowspan=2)
axSal.set_title('Salinity')

sectdistA, zA, salA, rhoA, yeardayA, latAS, lonAS = ReadAtlantisSection(7)

cmapsal = 'plasma'
conts = np.linspace(35.25, 36.75, 50)
im = axSal.contourf(sectdistA, zA, salA, conts, extend='both',cmap=cmapsal)
for c in im.collections:
    c.set_edgecolor("face")
cb = plt.colorbar(im, ax=axSal)
cb.set_ticks((conts[0], 36, conts[-1]))
cb.solids.set_edgecolor("face")

axSal.contour(sectdistA, zA, rhoA, 20, colors='k')
axSal.set_ylim((0, 200))
axSal.invert_yaxis()
axSal.set_ylabel('z [m]')
axSal.set_xlabel('Cross-stream distance [km]')

#Highlight Atlantis Section
axsst.plot(lonAS, latAS, linewidth=4, color='red')

plt.subplots_adjust(wspace=1, hspace=0.7)

#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/ObsOverview.pdf', bbox_inches='tight')

#%%
cmap = 'RdYlBu_r'
tliml = 12
tlimh = 23
plt.figure()
plt.pcolor(longrid[lonliml:lonlimr], latgrid[latliml:latlimr], np.transpose(sst[lonliml:lonlimr, latliml:latlimr]), cmap=cmap, vmin = tliml, vmax=tlimh)
plt.plot(flons, flats, linewidth=4)
plt.plot(lonK, latK, linewidth=4)
plt.plot(lonA, latA, linewidth=4)

indrelease = np.argmin(np.abs(fyds - 64.86))
plt.plot(flons[indrelease], flats[indrelease], marker='x')
plt.colorbar()
plt.xlim(-66, -62)
plt.ylim(38, 40)