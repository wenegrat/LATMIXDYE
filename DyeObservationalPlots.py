
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 15:24:57 2019

@author: jacob
"""

import h5py
import scipy.io as spio
import matplotlib.pyplot as plt
import numpy as np
from cmocean import cm as cmo
from scipy.interpolate import griddata
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#from matplotlib.mlab import griddata
import datetime as dt
import gsw
import scipy.interpolate as interp
from matplotlib.gridspec import GridSpec
#%%
plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 14
plt.rcParams['contour.negative_linestyle'] = 'solid'


#%% LOAD DYE DATA
#filename = '/home/jacob/dedalus/LATMIX/LatMix_2012/inj_data_2012'
#matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)
#matfile = spio.loadmat(filename,struct_as_record=True, squeeze_me=True)

filename = '/home/jacob/dedalus/LATMIX/LatMix_2012/my_triaxus_SI2'
matfile = spio.loadmat(filename,struct_as_record=True, squeeze_me=True)


depth = matfile['DEPTH']
fluorppb = matfile['FLUORPPB']
shiplog = matfile['SHIPLOG']
rho = matfile['PDENS']
jday = matfile['JDAY']
lat = matfile['LAT']
lon = matfile['LON']
II = matfile['II'] # This indexes individual gulf stream crossings

fluorppb[np.isnan(fluorppb)] = 0
fluorppb[fluorppb<0] = 0
fluorppb[fluorppb==0] = 1e-10
nd, ns = shiplog.shape

#%% LOAD TIMESERIES DATA

jday_ts = matfile['jday_ts']
fluorppb_ts = matfile['fluorPPB_ts']
depth_ts = matfile['depth_ts']
lat_ts = matfile['lat_ts']
lon_ts = matfile['lon_ts']
rho_ts = matfile['pdens_ts']
T_ts = matfile['T_ts']
S_ts = matfile['S_ts']
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

#sf = 10
#poly = np.polyfit(fyd,flat,6)
#flat_s = np.poly1d(poly)(fyd)
#
#window_size, poly_order = 51, 3
#flat_s = savgol_filter(flat, window_size, poly_order, deriv=0)
#flon_s = savgol_filter(flon, window_size, poly_order, deriv=0)
fvellat = np.gradient(flats)/np.gradient(fyds*86400)
fvellon = np.gradient(flons)/ np.gradient(fyds*86400)
fdir = np.arctan2(fvellat, fvellon)

filename = '/home/jacob/dedalus/LATMIX/FloatData/Env.mat'
matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)
fT = matfile['T']
fS = matfile['S']
fT = np.mean(fT, axis=-1)# Top and bottom sensors
fS = np.mean(fS, axis=-1)
fydr = matfile['yd'] - 1
fP = matfile['P']
#%% SCATTER IN T-S space
plt.figure()
norm = np.nanmax(fluorppb_ts)
mask = fluorppb_ts/norm > 1e-3
plt.scatter(S_ts[mask], T_ts[mask], c=np.log10(fluorppb_ts[mask]/norm), s=1)

#%% Calculated depth integrated 
# Might be easier to do from gridded data

# NB. Need to figure out how to calculate across-stream distance for each cast, otherwise
# Nb. Also need to figure out how to interpolate/plot surface density for each of these.

# 1 ) At each point, calculate distance from float
# 2 ) Using instantaneous float velocity calculate cross-stream distance

fluor_int = np.zeros((ns,))
lat_int = np.zeros((ns,))
lon_int = np.zeros((ns,))
jday_int = np.zeros((ns,))
surf_rho = np.zeros((ns,))
shiplog_int = np.zeros((ns,))
zcom_rho = np.zeros((ns,))
for i in range(0, ns):
    mask = np.isfinite(fluorppb[:,i])
    fluor_int[i] = integrate.trapz(fluorppb[mask,i], x=depth[mask,i])
    shiplog_int[i] = np.nanmean(shiplog[mask,i])
    lat_int[i] = np.nanmean(lat[:,i])
    lon_int[i] = np.nanmean(lon[:,i])
    jday_int[i] = np.nanmean(jday[:,i])
    ind = np.where(np.isfinite(rho[:,i]))[0][0]
    zcom = integrate.trapz(fluorppb[mask,i]*depth[mask,i], x=depth[mask,i])/fluor_int[i]
    surf_rho[i] = rho[ind,i]
    ind = np.argmin(np.abs(depth[:,i] - zcom))
    zcom_rho[i] = rho[ind,i]
# calc distance relative to float
iflat = np.interp(jday_int, fyds, flats)
iflon = np.interp(jday_int, fyds, flons)
distlat = (lat_int-iflat)*111e3*np.cos(lat_int)
distlon = (lon_int-iflon)*111e3

fdist = np.sign(lat_int-iflat)*np.sqrt(((lat_int - iflat)*111e3*np.cos(lat_int))**2 + ((lon_int - iflon)*111e3)**2)
#%% HISTO IN T-S Space 3 panel 

cmap = 'gnuplot'
cl = [-5, -1]
norm = np.max(fluorppb_ts)
mask = (fluorppb_ts/norm>10**(cl[0])) & (jday_ts<68)
#mask = jday_ts<90
xedges = np.linspace(34.5, 37, 150)
yedges = np.linspace(13, 22, 150)
X, Y = np.meshgrid(xedges, yedges)
SA = gsw.SA_from_SP(X, 0, -66, 39)
CT = gsw.CT_from_t(SA, Y, 0)
R = gsw.rho(SA, CT, 0)
R = R - 1000

SA = gsw.SA_from_SP(fS, fP, -66, 39)
CT = gsw.CT_from_t(SA, fT, fP)
Rf = gsw.rho(SA, CT, 0)
Rf = Rf - 1000
norm = np.sum(fluorppb_ts[mask])
norm = 1
H, xedges, yedges = np.histogram2d(S_ts[mask], T_ts[mask], weights=(fluorppb_ts[mask]/norm),   bins=(xedges, yedges), density=False)
fig = plt.figure(figsize=(12.78, 9.41))

# TIME-PLOT
axhistot = plt.subplot2grid((3,4), (0,2), rowspan = 2, colspan=2)

Htime, xedges, yedges = np.histogram2d(S_ts[mask], T_ts[mask], weights=(fluorppb_ts[mask]*jday_ts[mask]),   bins=(xedges, yedges), density=False)
Hcount, xedges, yedges = np.histogram2d(S_ts[mask], T_ts[mask],   bins=(xedges, yedges), density=False)
Hf, xedges, yedges = np.histogram2d(S_ts[mask], T_ts[mask],  weights=(fluorppb_ts[mask]), bins=(xedges, yedges), density=False)

Havgtime = Htime/(Hf) # This is the concentration weighted average time
#Havgtime = Htime/Hcount
cmapt = 'gnuplot'
ix = axhistot.pcolor( xedges, yedges, np.transpose(Havgtime),  cmap=cmapt, vmin=64.87, vmax=68)
ix.set_edgecolor('face')
#axhisto.grid()
conts = np.linspace(20, 30, 21)
CL = axhistot.contour(X, Y, R,conts, colors='k', linestyles='dashed', linewidths=0.5)
axhistot.clabel(CL, inline=1, fontsize=10, fmt='%1.1f')
axhistot.set_xlabel('Salinity [psu]')
axhistot.set_ylabel('Temperature [$^\circ$ C]')
si = np.argmin(np.abs(fydr-65))
ei = np.argmin(np.abs(fydr-68))
dstep = (fydr[-1] - fydr[0])/fydr.size
span = int(1/6/dstep)
lims = range(si, ei, span)
axhistot.scatter(fS[lims], fT[lims], c = fydr[lims], marker='d', edgecolor='w', s=40, cmap=cmapt, vmin=64.87, vmax=68)
axhistot.text(34.55, 20.85 ,'b)', color='k', size=20, bbox=dict(facecolor='w', edgecolor='k'))

cbaxes = inset_axes(axhistot, width="40%", height="1.5%", loc=4, borderpad=1.25) 

cb = plt.colorbar(ix,cax=cbaxes,orientation='horizontal',)
cb.ax.xaxis.set_ticks_position('top')
cb.set_ticks([65, 66, 67, 68])
#cb = plt.colorbar(ix, extend='both')
cb.set_label('Yearday')

# COLORED HISTOGRAM
axhisto = plt.subplot2grid((3,4), (0,0), rowspan = 2, colspan=2)
ix = axhisto.pcolor( xedges, yedges, np.transpose(np.log10(H/np.sum(H))), vmin=cl[0], vmax=cl[1], cmap=cmap)
ix.set_edgecolor('face')
conts = np.linspace(20, 30, 21)
CL = axhisto.contour(X, Y, R,conts, colors='k', linestyles='dashed', linewidths=0.5)

axhisto.clabel(CL, inline=1, fontsize=10, fmt='%1.1f')
axhisto.set_xlabel('Salinity [psu]')
axhisto.set_ylabel('Temperature [$^\circ$ C]')
#si = np.argmin(np.abs(fydr-65))
#dstep = (fydr[-1] - fydr[0])/fydr.size
#span = int(1/6/dstep)
#lims = range(si, fydr.size, span)
axhisto.scatter(fS[lims], fT[lims], marker='d', color='cornflowerblue', edgecolor='w', s=40)
axhisto.text(34.55, 20.85 ,'a)', color='k', size=20, bbox=dict(facecolor='w', edgecolor='k'))

#cbaxes = inset_axes(axhisto, width="1%", height="60%", loc=6) 
#cb = plt.colorbar(ix, cax=cbaxes, extend='both')
cbaxes = inset_axes(axhisto, width="40%", height="1.5%", loc=4, borderpad=1.25) 

cb = plt.colorbar(ix,cax=cbaxes,orientation='horizontal', extend='both')
cb.ax.xaxis.set_ticks_position('top')
cb.set_ticks([-5,  -3,  -1])
#cb = plt.colorbar(ix, extend='both')
cb.set_label('log$_{10}$(Relative Frequency)')


## SURFACE PLOT
cl=[-3, 0]
axsurf = plt.subplot2grid((3,4), (2,0), rowspan=1, colspan=4)
normscatter = np.nanmax(fluor_int)
axsurf.grid()
mask = fluor_int/normscatter>1e-5
ix = axsurf.scatter(jday_int[mask], zcom_rho[mask], c=np.log10(fluor_int[mask]/normscatter), cmap=cmap, vmin=cl[0], vmax=cl[1])
axsurf.scatter(fydr[lims], Rf[lims], marker='d', color='cornflowerblue', edgecolor='w', s=60)
axsurf.set_ylabel('Density, $\sigma_0$ [kg m$^{-3}$]')
axsurf.set_xlabel('Yearday')
axsurf.set_xlim(64.5, 68)
axsurf.set_ylim(25, 27)
axsurf.text(64.55, 26.675 ,'c)', color='k', size=20, bbox=dict(facecolor='w', edgecolor='k'))

cbaxes = inset_axes(axsurf, width="40%", height="3%", loc=4, borderpad=1.1) 

cb = plt.colorbar(ix,cax=cbaxes,orientation='horizontal', extend='min')
cb.set_label('log$_{10}$(Normalized Concentration)', fontsize=12)
cb.set_ticks([-3, -2, -1, 0])
cb.ax.tick_params(labelsize=12, pad=-1)
#cb.ax.xaxis.set_label_position('top')
cb.ax.xaxis.set_ticks_position('top')
plt.subplots_adjust(wspace=0.4, hspace=0.4)

#axs['5'].contour(jd, sd, grid_surf_rho, 5)
#axs['5'].set_xlim([64.75, 67])
#axsurf.set_ylim(-10, 10)
#plt.scatter(S_ts[mask], T_ts[mask], c=np.log10(fluorppb_ts[mask]/norm), s=1)

#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/DyeDensity.pdf', bbox_inches='tight')




#%% HISTO IN T-S Space

cmap = 'gnuplot'
cl = [-3, 0]

mask = (fluorppb_ts/norm>1e-5) & (jday_ts<68)
#mask = jday_ts<90
xedges = np.linspace(34.5, 37, 100)
yedges = np.linspace(13, 22, 100)
X, Y = np.meshgrid(xedges, yedges)
SA = gsw.SA_from_SP(X, 0, -66, 39)
CT = gsw.CT_from_t(SA, Y, 0)
R = gsw.rho(SA, CT, 0)
R = R - 1000

SA = gsw.SA_from_SP(fS, fP, -66, 39)
CT = gsw.CT_from_t(SA, fT, fP)
Rf = gsw.rho(SA, CT, 0)
Rf = Rf - 1000
norm = np.sum(fluorppb_ts[mask])
norm = 1
H, xedges, yedges = np.histogram2d(S_ts[mask], T_ts[mask], weights=(fluorppb_ts[mask]/norm),   bins=(xedges, yedges), density=False)
fig = plt.figure(figsize=(7.26, 9.4))
axhisto = plt.subplot2grid((3,2), (0,0), rowspan = 2, colspan=2)

ix = axhisto.pcolor( xedges, yedges, np.transpose(np.log10(H/np.sum(H))), vmin=-6, vmax=-1, cmap=cmap)
ix.set_edgecolor('face')
#axhisto.grid()
conts = np.linspace(20, 30, 21)
CL = axhisto.contour(X, Y, R,conts, colors='k', linestyles='dashed', linewidths=0.5)
axhisto.clabel(CL, inline=1, fontsize=10, fmt='%1.1f')
axhisto.set_xlabel('Salinity [psu]')
axhisto.set_ylabel('Temperature [$^\circ$ C]')
si = np.argmin(np.abs(fydr-65))
dstep = (fydr[-1] - fydr[0])/fydr.size
span = int(1/6/dstep)
lims = range(si, fydr.size, span)
axhisto.scatter(fS[lims], fT[lims], marker='d', color='cornflowerblue', edgecolor='w', s=40)

#cbaxes = inset_axes(axhisto, width="1%", height="60%", loc=6) 
#cb = plt.colorbar(ix, cax=cbaxes, extend='both')
cbaxes = inset_axes(axhisto, width="40%", height="1.5%", loc=4, borderpad=1.25) 

cb = plt.colorbar(ix,cax=cbaxes,orientation='horizontal', extend='both')
cb.ax.xaxis.set_ticks_position('top')
#cb = plt.colorbar(ix, extend='both')
cb.set_label('log$_{10}$(Relative Frequency)')
axsurf = plt.subplot2grid((3,2), (2,0), rowspan=1, colspan=2)
normscatter = np.nanmax(fluor_int)
axsurf.grid()
mask = fluor_int/normscatter>1e-5
ix = axsurf.scatter(jday_int[mask], zcom_rho[mask], c=np.log10(fluor_int[mask]/normscatter), cmap=cmap, vmin=cl[0], vmax=cl[1])
axsurf.scatter(fydr[lims], Rf[lims], marker='d', color='cornflowerblue', edgecolor='w', s=60)
axsurf.set_ylabel('Density, $\sigma_0$ [kg m$^{-3}$]')
axsurf.set_xlabel('Yearday')
axsurf.set_xlim(64.5, 68)
axsurf.set_ylim(25., 27)
cbaxes = inset_axes(axsurf, width="40%", height="3%", loc=4, borderpad=1.1) 

cb = plt.colorbar(ix,cax=cbaxes,orientation='horizontal', extend='min')
cb.set_label('log$_{10}$(Normalized Concentration)', fontsize=12)
cb.set_ticks([-3, -2, -1, 0])
cb.ax.tick_params(labelsize=12, pad=-1)
#cb.ax.xaxis.set_label_position('top')
cb.ax.xaxis.set_ticks_position('top')
plt.subplots_adjust(wspace=0.2, hspace=0.4)

#axs['5'].contour(jd, sd, grid_surf_rho, 5)
#axs['5'].set_xlim([64.75, 67])
#axsurf.set_ylim(-10, 10)
#plt.scatter(S_ts[mask], T_ts[mask], c=np.log10(fluorppb_ts[mask]/norm), s=1)

#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/DyeDensity.pdf', bbox_inches='tight')




#%%
plt.figure()
ll = np.argmin(np.abs(jday_ts - 66))
plt.scatter(jday_ts[0:ll], rho_ts[0:ll], s=[], c=np.log10(fluorppb_ts[0:ll]))

#%%
plt.figure()
ll = np.argmin(np.abs(jday_int - 66))
plt.scatter(jday_int[0:ll], surf_rho[0:ll], s=[], c=np.log10(fluor_int[0:ll]))
#%% INTEGRATED BY CROSSING

nc, nl = II.shape
survdist = np.zeros((ns,))
weightedrho = np.zeros((nc,))
ydc = np.zeros((nc,))
shiplogc = np.zeros((nd, ns))
for i in range(0, nc-1):
    span = range(II[i,0], II[i,1])
#    distspan = fdist[span]
#    distind = np.argmin(np.abs(distspan))
    
    meanfloatlat = np.interp(np.mean(jday_int[span]), fyds, flats)
    meanfloatlon = np.interp(np.mean(jday_int[span]), fyds, flons)
    distlat = (lat_int[span]-meanfloatlat)*111e3*np.cos(lat_int[span])
    distlon = (lon_int[span]-meanfloatlon)*111e3
    distspan = np.sign(distlat)*np.sqrt(distlat**2 + distlon**2)
    distind = np.argmin(np.abs(distspan))

    tempdist = distspan - distspan[distind] # This is distance from the closest point of approach to the float for a survey
    survang = np.arctan2((lat_int[span[-1]] - lat_int[span[0]]), (lon_int[span[-1]] - lon_int[span[0]]))
    
    fvellatm = np.mean(np.interp(jday_int[span], fyds, fvellat))
    fvellonm = np.mean(np.interp(jday_int[span], fyds, fvellon))
    fang = np.arctan2(fvellatm, fvellonm)
#    print(survang*180/np.pi)
    angdif = survang - fang
#    angdif = (angdif+np.pi) % (2*np.pi) - np.pi
    print(angdif*180/np.pi)
    factor = np.sin(angdif)
    
    survdist[span] = (shiplog_int[span] - shiplog_int[span][distind])*factor
#    survdist[span] = tempdist*np.abs(factor)
    shiplogc[:,span] = (shiplog[:,span] - shiplog[:,span][:,distind][:,np.newaxis])*factor
#    print(survang*180/np.pi)
    ftemp = fluorppb[:,span]
    rtemp = rho[:,span]
    ftemp[~np.isfinite(ftemp)] = 0
    rtemp[~np.isfinite(rtemp)] = 0
    sltemp = shiplogc[:,span]
    sltemp[~np.isfinite(sltemp)] = 0
    weightedrho[i] = np.trapz(np.trapz(ftemp*rtemp, x=sltemp, axis=-1), x=depth[:,span[0]], axis=0)/np.trapz(np.trapz(ftemp, x=sltemp, axis=-1), x=depth[:,span[0]], axis=0)
    ydc[i] = np.nanmean(jday[:,span])
plt.figure()
plt.scatter(lon_int[span],lat_int[span], c=fluor_int[span])
plt.scatter(lon_int[span],lat_int[span], c=distspan)

plt.colorbar()
plt.scatter(np.interp(np.mean(jday_int[span]), fyds, flons),np.interp(np.mean(jday_int[span]), fyds, flats), marker='x')
plt.scatter(lon_int[span][distind], lat_int[span][distind], c='r')
plt.quiver(np.interp(np.mean(jday_int[span]), fyds, flons),np.interp(np.mean(jday_int[span]), fyds, flats), fvellonm, fvellatm)
plt.axes().set_aspect('equal', 'datalim')

plt.figure()
plt.plot(survdist[span], fluor_int[span], marker='x')

#%% INTERPOLATE SURFACE DENSITY
jd = np.linspace(jday_int[0], jday_int[-1], 100)
sd = np.linspace(-10, 10, 100)

gjd, gsd = np.meshgrid(jd, sd)
surf_rho[surf_rho < 24] = np.nan
mask = np.isfinite(jday_int + survdist + surf_rho)
grid_surf_rho = griddata((jday_int[mask], survdist[mask]), surf_rho[mask], (gjd, gsd), method='nearest')


#%% Plot 4 panel
##################### NOT INTERPOLATED #####################################
cmap = 'gnuplot'
conts = np.linspace(-3, 0, 13)

crossings = [6, 8, 12, 14] # pcik the sections 
crossings = [4, 6, 8, 10]
crossings = [4,  6, 8, 12]
#crossings = [0, 1, 2, 3]
xcrosslims = [[-5, 5], [-5, 5], [-4, 6], [-7, 3]]
#crossings = [4, 6, 7, 8]
norm = np.nanmax(fluorppb)


fig=plt.figure(figsize=(12.2,6.6))

gs=GridSpec(2,2) # 3 rows, 2 columns
axs = {}
axs['1']=fig.add_subplot(gs[0,0]) # First row, first column
axs['2']=fig.add_subplot(gs[0,1]) # First row, second column
axs['3']=fig.add_subplot(gs[1,0]) # First row, first column
axs['4']=fig.add_subplot(gs[1,1]) # First row, second column

cl = [0.001, 1]
cl = [-3, 0]
rhoc = np.linspace(20, 30, 101)+0.06
dyeconts = np.linspace(-3, 0, 13)
# First section
for i in range(1, len(crossings)+1):
    
    span = range(II[crossings[i-1],0], II[crossings[i-1],1])
    X = shiplog[:,span] - shiplog[:,span[0]][:,np.newaxis]
    X = shiplogc[:,span]
    Y = depth[:,span]
    C = fluorppb[:,span]/norm
    R = rho[:,span]
    im = axs[str(i)].pcolor(X, Y, np.log10(C), vmin=cl[0], vmax=cl[1], cmap=cmap)
#    im = axs[str(i)].contourf(X, Y, np.log10(C), dyeconts, vmin=cl[0], vmax=cl[1], cmap=cmap)
    im.set_edgecolor('face')
    axs[str(i)].contour(X,Y, R,rhoc, colors='w')
    axs[str(i)].set_ylim(0, 150)
    axs[str(i)].invert_yaxis()
    axs[str(i)].set_xlabel('km')
    axs[str(i)].set_ylabel('Pressure [db]')
    axs[str(i)].set_title('Yearday: %2.2f' %np.nanmean(jday[:,span]))
    axs[str(i)].set_xlim(xcrosslims[i-1])
#    axs[str(i)].grid()


#test float data
#axs['5'].plot(fyd, flat)
# ADD COLORBAR
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.91, 0.25, 0.01, 0.5])
cb = fig.colorbar(im, cax=cbar_ax, extend='min')
cb.set_ticks([-3, -2, -1, 0])
cb.set_label('log$_{10}$(Normalized Concentration)')
cb.solids.set_edgecolor("face")
plt.subplots_adjust(wspace=0.2, hspace=0.4)

#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/DyeObsOverview.pdf', bbox_inches='tight')

#%% Plot with all surveys
##################### NOT INTERPOLATED #####################################
cmap = 'gnuplot'
conts = np.linspace(-3, 0, 13)

crossings = [6, 8, 12, 14] # pcik the sections 
crossings = [4, 6, 8, 10]
crossings = [4,  6, 8, 12]
crossings = [0, 1, 2, 3]
xcrosslims = [[-5, 5], [-5, 5], [-4, 6], [-7, 3]]
#crossings = [4, 6, 7, 8]
norm = np.nanmax(fluorppb)


fig=plt.figure(figsize=(12.5,10))

gs=GridSpec(3,2) # 3 rows, 2 columns
axs = {}
axs['1']=fig.add_subplot(gs[0,0]) # First row, first column
axs['2']=fig.add_subplot(gs[0,1]) # First row, second column
axs['3']=fig.add_subplot(gs[1,0]) # First row, first column
axs['4']=fig.add_subplot(gs[1,1]) # First row, second column
axs['5']=fig.add_subplot(gs[2,:]) # Second row, span all columns

cl = [0.001, 1]
cl = [-3, 0]
rhoc = np.linspace(20, 30, 101)+0.06
dyeconts = np.linspace(-3, 0, 13)
# First section
for i in range(1, len(crossings)+1):
    
    span = range(II[crossings[i-1],0], II[crossings[i-1],1])
    X = shiplog[:,span] - shiplog[:,span[0]][:,np.newaxis]
    X = shiplogc[:,span]
    Y = depth[:,span]
    C = fluorppb[:,span]/norm
    R = rho[:,span]
    im = axs[str(i)].pcolor(X, Y, np.log10(C), vmin=cl[0], vmax=cl[1], cmap=cmap)
#    im = axs[str(i)].contourf(X, Y, np.log10(C), dyeconts, vmin=cl[0], vmax=cl[1], cmap=cmap)
    im.set_edgecolor('face')
    axs[str(i)].contour(X,Y, R,rhoc, colors='w')
    axs[str(i)].set_ylim(0, 150)
    axs[str(i)].invert_yaxis()
    axs[str(i)].set_xlabel('km')
    axs[str(i)].set_ylabel('Pressure [db]')
    axs[str(i)].set_title('Yearday: %2.2f' %np.nanmean(jday[:,span]))
    axs[str(i)].set_xlim(xcrosslims[i-1])
#    axs[str(i)].grid()

# ADD INTEGRATED PLOT
normscatter = np.nanmax(fluor_int)
axs['5'].grid()
axs['5'].scatter(jday_int, survdist, c=np.log10(fluor_int/normscatter), cmap=cmap, vmin=cl[0], vmax=cl[1])
axs['5'].set_ylabel('Cross-front distance [km]')
axs['5'].set_xlabel('Yearday')

#axs['5'].contour(jd, sd, grid_surf_rho, 5)
#axs['5'].set_xlim([64.75, 67])
axs['5'].set_ylim(-10, 10)

#test float data
#axs['5'].plot(fyd, flat)
# ADD COLORBAR
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.91, 0.25, 0.01, 0.5])
cb = fig.colorbar(im, cax=cbar_ax, extend='min')
cb.set_ticks([-3, -2, -1, 0])
cb.set_label('log$_{10}$(Normalized Concentration)')
cb.solids.set_edgecolor("face")
plt.subplots_adjust(wspace=0.2, hspace=0.4)

#plt.savefig('/home/jacob/Dropbox/WeeklyReports/190801/DyeObsOverview.pdf', bbox_inches='tight')
#%%
plt.figure()
plt.scatter(jday_int, lat_int, c=np.log10(fluor_int/normscatter), cmap=cmap, vmin=cl[0], vmax=cl[1])
plt.plot(fyds, flats)
#%% Make plot of surface dye density and dye concentration
################### INTERPOLATED (SLOW) ####################################
cmap = 'gnuplot'
conts = np.linspace(-3, 0, 13)

crossnum = 12 # yearday 65.52
#crossnum = 6 # first day there is obs 65.24


norm = np.nanmax(fluorppb)


span = range(II[crossnum,0], II[crossnum,1])
X = shiplog[:,span]
Y = depth[:,span]
C = fluorppb[:,span]
# indices where nobody is nan
inds = np.logical_not(np.isnan(X) | np.isnan(Y) | np.isnan(C))
X_notnan = X[inds]
Y_notnan = Y[inds]
C_notnan = C[inds]

# construct new mesh
X_vals = np.unique(X[inds])
Y_vals = np.unique(Y[inds])
X_plot,Y_plot = np.meshgrid(X_vals,Y_vals)

# use nearest-neighbour interpolation to match the two meshes
C_plot = interp.griddata(np.array([X_notnan,Y_notnan]).T,C_notnan,(X_plot,Y_plot),method='nearest')
# fill in the nans in C
C_plot[np.logical_not(inds)] = np.nan
#%%
fig, ax = plt.subplots()
n_levels = 100  # number of contour levels to plot
ax.contourf(X_plot,Y_plot,C_plot,n_levels)
#%%

plt.figure()

#plt.pcolor(jday, shiplog, rho)
#plt.contourf( shiplog[:,span],-depth[:,span], (np.log10(fluorppb[:,span]/norm)),conts, extend='both',cmap=cmap)
#plt.contourf(X_plot, -Y_plot, np.log10(C_plot/norm),conts, extend='both',cmap=cmap)
plt.pcolor(X_plot, -Y_plot, np.log10(C_plot/norm),cmap=cmap)

#plt.clim(0, 3)
plt.colorbar()
plt.title(np.nanmean(jday[:,span]))
#plt.contour( shiplog[:,span],-depth[:,span], rho[:,span], 10,colors='w')

