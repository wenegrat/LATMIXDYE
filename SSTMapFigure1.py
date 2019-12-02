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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from ReadAtlantis_ALL import ReadAtlantis_ALL
from ReadAtlantis import ReadAtlantisSection
from ReadKnorr_ALL import ReadKnorr_ALL

import cartopy.crs as ccrs
import cartopy.feature as cfeature
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


depth = matfile['DEPTH']
fluorppb = matfile['FLUORPPB']
shiplog = matfile['SHIPLOG']
rho = matfile['PDENS']
jday = matfile['JDAY']
lat = matfile['LAT']
lon = matfile['LON']
II = matfile['II'] # This indexes individual gulf stream crossings
crossings = [4,  6, 8, 12]

# Get the sections used in the dye observatoinal plot
span = range(II[crossings[0],0], II[crossings[0],1])
latl1 = np.nanmean(lat[:,span],axis=0)
lonl1 = np.nanmean(lon[:,span], axis=0)

span = range(II[crossings[1],0], II[crossings[1],1])
latl2 = np.nanmean(lat[:,span],axis=0)
lonl2 = np.nanmean(lon[:,span], axis=0)

span = range(II[crossings[2],0], II[crossings[2],1])
latl3 = np.nanmean(lat[:,span],axis=0)
lonl3 = np.nanmean(lon[:,span], axis=0)

span = range(II[crossings[3],0], II[crossings[3],1])
latl4 = np.nanmean(lat[:,span],axis=0)
lonl4 = np.nanmean(lon[:,span], axis=0)
#%% LOAD ATLANTIS TRACK
latA, lonA, timeA = ReadAtlantis_ALL()
sectdistA, zA, salA, rhoA, yeardayA, latAS, lonAS = ReadAtlantisSection(7)
zA  = -zA
#%%
cmap = 'RdYlBu_r'
tliml = 13
tlimh = 23

cols = 7

fig = plt.figure(figsize=(12.22, 6.16))

# SST PLOT

aspratio = 1/np.cos(39*np.pi/180) # > 1 indicates x coordinate is smaller than y coordinate

axsst = plt.subplot2grid((3,cols), (0,0), rowspan = 2, colspan=4)
isst = axsst.pcolor(longrid[lonliml-1:lonlimr+1], latgrid[latliml-1:latlimr+1], np.transpose(sst[lonliml-1:lonlimr+1, latliml-1:latlimr+1]), cmap=cmap, vmin = tliml, vmax=tlimh)
isst.set_edgecolor('face')
axsst.plot(flons, flats, linewidth=2,  color='k', label='Float')
timeAind = np.argmin(np.abs(timeA - 64.86))
axsst.plot(lonA[timeAind:], latA[timeAind:], linewidth=2, color = '#2ca02c', label='R/V Atlantis')
#axsst.plot(lonAS[0], latAS[0], linewidth=1, color='#2ca02c', marker='o')  #Highlight Atlantis Section
#axsst.plot(lonAS[-1], latAS[-1], linewidth=1, color='#2ca02c', marker='o')  #Highlight Atlantis Section
l1, = axsst.plot(lonAS, latAS, color='w', linestyle='dashed')
l1.set_dashes([1, 1])
           
axsst.plot(lonK, latK, linewidth=2, color = '#1f77b4', label='R/V Knorr')
# mark Knorr sections used
l1, = axsst.plot(lonl1, latl1, color = 'w', linestyle='dashed')
l1.set_dashes([1, 1])

l1, = axsst.plot(lonl2, latl2, color = 'w', linestyle='dashed')
l1.set_dashes([1, 1])
    
l1, = axsst.plot(lonl3, latl3, color = 'w', linestyle='dashed')
l1.set_dashes([1, 1])

l1, = axsst.plot(lonl4, latl4, color = 'w', linestyle='dashed')
l1.set_dashes([1, 1])
                                        
indrelease = np.argmin(np.abs(fyds - 64.86))
axsst.plot(flons[indrelease], flats[indrelease], marker='o', color='w')
axsst.annotate('Dye release', xytext=(-65.975, 38.75), xy=(flons[indrelease], flats[indrelease]), color='w',
            arrowprops=dict(arrowstyle='->', color='w'),
            )
cbaxes = inset_axes(axsst, width="40%", height="3%", loc=9) 
cb = plt.colorbar(isst, cax=cbaxes, orientation='horizontal', extend='both')
cb.set_label('SST [$^{\circ}$ C]', color='w', fontsize=13)
cb.outline.set_edgecolor('w')
cb.set_ticks([tliml, (tlimh-tliml)/2+tliml, tlimh])
cb.ax.xaxis.set_tick_params(color='w')
cb.ax.tick_params(labelsize=12)
plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='w') # set colorbar  
axsst.text( -65.93, 39.8 ,'a)', color='k', size=20, bbox=dict(facecolor='w', edgecolor='k'))
#plt.colorbar(isst, cax = axsst)
axsst.set_xlim(-66, -62)
axsst.set_ylim(38, 40)
axsst.set_aspect(aspratio)
axsst.set_ylabel('Latitude')
axsst.set_xlabel('Longitude')
axsst.grid(linestyle='dashed')
axsst.legend(loc=4, ncol=3, fontsize=11, handletextpad=0.2, columnspacing=0.8)


# SURFACE FORCING PLOT
axflux = plt.subplot2grid((3,cols), (2,1), colspan=cols-2)
color = 'tab:blue'
axflux.plot(yearday, tmag, color=color, linewidth=2)
axflux.tick_params(axis='y', labelcolor=color)
axflux.set_ylabel('Wind-stress\n magnitude\n [N m$^{-2}$]', color=color)  # we already handled the x-label with ax1
axflux.set_ylim((0, 1.2))
axflux.set_yticks([0, 0.3, 0.6, 0.9,  1.2])
axflux.set_xlim((64.5, 66))
axflux.set_xlabel('yearday')

axflux2 = axflux.twinx()
axflux.axvline(64.86, linestyle='dashed', color='k')
color = 'tab:red'
axflux2.plot(yearday, qnet, color=color, linewidth=2)
axflux2.set_ylabel('Net surface\n heat flux\n [W m$^{-2}$]', color=color)  # we already handled the x-label with ax1
axflux2.tick_params(axis='y', labelcolor=color)
axflux2.set_ylim((0, 1200))
axflux2.set_yticks([0, 300, 600, 900,  1200])
axflux.grid(linestyle='dashed')
axflux.text( -0.22, 0.75 ,'c)', color='k', size=20, bbox=dict(facecolor='w', edgecolor='k'),transform=axflux.transAxes )

# SALINITY SECTION
axSal = plt.subplot2grid((3,cols), (0,4), colspan=cols-4, rowspan=2)

cmapsal = cmo.haline
#cmapsal = 'RdYlBu_r'
conts = np.linspace(35.25, 36.75, 50)
im = axSal.contourf(sectdistA, zA, salA, conts, extend='both',cmap=cmapsal)
for c in im.collections:
    c.set_edgecolor("face")
cb = plt.colorbar(im, ax=axSal)
cb.set_ticks((conts[0], 36, conts[-1]))
cb.solids.set_edgecolor("face")
cb.set_label('Salinity [psu]')

axSal.contour(sectdistA, zA, rhoA, 20, colors='k')
axSal.set_ylim((-200, 0))
axSal.set_ylabel('z [m]')
axSal.set_xlabel('Cross-stream distance [km]')

axSal.text( -4.5, -20 ,'b)', color='k', size=20, bbox=dict(facecolor='w', edgecolor='k'))

plt.subplots_adjust(wspace=1, hspace=0.7)

#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/ObsOverview.pdf', bbox_inches='tight')

#%%

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(-45, 45))

ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND, edgecolor='black')
#
ax.set_global()
#ax.gridlines()


#ax.pcolormesh(longrid[:], latgrid[:], np.transpose(sst[:,:]), cmap=cmap, vmin = tliml, vmax=tlimh, transform = ccrs.PlateCarree())

tf = ccrs.Geodetic()
ax.plot([-66, -62], [38, 38],
         color='r', linestyle='-',
         transform=tf,
         )
ax.plot([-66, -62], [40, 40],
         color='r', linestyle='-',
         transform=tf,
         )
ax.plot([-66, -66], [38, 40],
         color='r', linestyle='-',
         transform=tf,
         )
ax.plot([-62, -62], [38, 40],
         color='r', linestyle='-',
         transform=tf,
         )
