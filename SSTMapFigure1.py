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