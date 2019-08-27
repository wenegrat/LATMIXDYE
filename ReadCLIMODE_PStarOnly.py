#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:28:10 2019

@author: jacob
"""

import h5py
import scipy.io as spio
import matplotlib.pyplot as plt
import numpy as np
from cmocean import cm as cmo
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
#from matplotlib.mlab import griddata
import gsw

plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 14
plt.rcParams['contour.negative_linestyle'] = 'solid'

from ReadAtlantis import ReadAtlantisSection

#%%
class CTDcast:
  def __init__(self):
    self.lat = None
    self.lon = None
    self.press = np.array([])
    self.theta = np.array([])
    self.sal = np.array([])
    self.oxygen = np.array([])
    self.silicate = np.array([])
    self.nitrate = np.array([])
    self.nitrite = np.array([])
    self.phosphate = np.array([])
    self.flags = None
  def setLat(self, lat, latm):
      self.lat = lat + latm/60
  def setLon(self, lon, lonm):
      self.lon = lon + lonm/60
  def getRho(self):
      SA = gsw.SA_from_SP(self.sal,self.press, self.lon, self.lat)
      CT = gsw.CT_from_t(SA, self.theta, self.press) # Is this potential temperature
      rho = gsw.rho(SA, CT, self.press)
      return rho
  def getRhoSorted(self):
      rhoraw = self.getRho()
      mask = np.isfinite(rhoraw)
      temp = rhoraw.copy()
      temp[mask] = np.sort(rhoraw[mask], axis=0)
      pressort = self.press.copy()
      pressort[mask] = np.sort(self.press[mask], axis=0)
      return temp, pressort
  def getDistance(self, startlat, startlon):
      latdeg = self.lat*np.pi/180
      londeg = self.lon*np.pi/180
      startlatdeg = startlat*np.pi/180
      startlondeg = startlon*np.pi/180
      
      d=2*np.arcsin(np.sqrt((np.sin((startlatdeg-latdeg)/2))**2 + 
                 np.cos(startlatdeg)*np.cos(latdeg)*(np.sin((startlondeg-londeg)/2))**2))
      dkm = 2*np.pi*6356.752*d/(2*np.pi)*np.sign(latdeg-startlatdeg)
      return dkm
#  def checkQualityFlags(self):
#      parseflags = list(self.flags)
#      if any(int(x)>2 for x in parseflags):
#          self.theta = -9 + 0*self.theta
#          self.sal = -9 + 0*self.theta
#          self.oxygen = -9 + 0*self.oxygen
#          self.silicate = -9 + 0*self.oxygen
#          self.nitrate = -9 + 0*self.nitrate
#          self.nitrite = -9 + 0*self.nitrite
#          self.phosphate = -9 + 0*self.phosphate
      
#%% LOAD THE DATA
f = open('/home/jacob/dedalus/LATMIX/CLIMODEDATA/33AT013.sum.txt', 'r')
fd = open('/home/jacob/dedalus/LATMIX/CLIMODEDATA/33AT013.sea.txt', 'r')
#data = f.read()
#f.close()

c1 = 0
stations = dict()

for line in f:
    c1 += 1
    line.strip()
    columns = line.split()
    if c1>4:
        lat = float(columns[8])
        latm = float(columns[9])
        lon = -float(columns[11])
        lonm = -float(columns[12])
        stnnum = columns[2]
        stations[stnnum] = CTDcast()
        stations[stnnum].setLat(lat, latm)
        stations[stnnum].setLon(lon, lonm)
#    print(repr(line))
    
c2 = 0
for line in fd:
    c2 +=1
    if c2>4:
        columns = line.split()
        stnnum = columns[0]
        if True: #int(columns[1]) == 1: #ignore second cast
            flags = list(columns[15])
            
            if all(int(x)<3 for x in flags):
                stations[stnnum].press   = np.concatenate((stations[stnnum].press, [float(columns[4])] ), axis=0)
                stations[stnnum].theta   = np.concatenate((stations[stnnum].theta, [float(columns[8])] ), axis=0)
                stations[stnnum].sal     = np.concatenate((stations[stnnum].sal, [float(columns[9])] ), axis=0)
                stations[stnnum].oxygen  = np.concatenate((stations[stnnum].oxygen, [float(columns[10])] ), axis=0)
                stations[stnnum].silicate= np.concatenate((stations[stnnum].silicate, [float(columns[11])] ), axis=0)
                stations[stnnum].nitrate = np.concatenate((stations[stnnum].nitrate, [float(columns[12])] ), axis=0)
                stations[stnnum].nitrite = np.concatenate((stations[stnnum].nitrite, [float(columns[13])] ), axis=0)
                stations[stnnum].phosphate= np.concatenate((stations[stnnum].phosphate, [float(columns[14])] ), axis=0)
            else:
                print('Bad Value, station num: ' + stnnum)
                badval = -9.0
                stations[stnnum].press   = np.concatenate((stations[stnnum].press, [badval] ), axis=0)
                stations[stnnum].theta   = np.concatenate((stations[stnnum].theta, [badval] ), axis=0)
                stations[stnnum].sal     = np.concatenate((stations[stnnum].sal, [badval] ), axis=0)
                stations[stnnum].oxygen  = np.concatenate((stations[stnnum].oxygen, [badval] ), axis=0)
                stations[stnnum].silicate= np.concatenate((stations[stnnum].silicate, [badval] ), axis=0)
                stations[stnnum].nitrate = np.concatenate((stations[stnnum].nitrate, [badval] ), axis=0)
                stations[stnnum].nitrite = np.concatenate((stations[stnnum].nitrite, [badval] ), axis=0)
                stations[stnnum].phosphate= np.concatenate((stations[stnnum].phosphate, [badval] ), axis=0)
sections = dict()
sections['1'] = [1,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
sections['2'] = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
sections['3'] = [25, 26, 27, 28, 29, 30, 31, 32, 33]
sections['4'] = [35, 36, 37, 38, 39, 40, 41, 42, 43]
sections['5'] = [45, 46, 47, 48, 49, 50]
#LOAD ADCP            
filename = '/home/jacob/dedalus/LATMIX/CLIMODEDATA/at13vmadcp.mat'
fadcp = spio.loadmat(filename)

#%%
# List all groups
print("Keys: %s" % fadcp.keys())
keys = list(fadcp.keys())[:]

adcpdepths = fadcp['depth'][0,:]
time = fadcp['txy'][:,0]
lat = fadcp['txy'][:,2]
lon = fadcp['txy'][:,1]

u = fadcp['u']
v = fadcp['v']
vtemp = v.copy()
vtemp[~np.isfinite(vtemp)] = 0
utemp = u.copy()
utemp[~np.isfinite(utemp)] = 0
ubar = np.mean(utemp, axis=-1)
vbar = np.mean(vtemp, axis=-1)
angle = np.arctan2(vbar, ubar)
velmag = np.sqrt(ubar**2 + vbar**2)


#%% GRID THE DATA
def loadSection(secnum):
    # First pick the stations:
    lats = []
    lons = []
    thetas= []
    press = []
    nitrates=[]
    phosphates = []
    rho = []
    oxygens = []
    distance = []
    
#    secnum = '1'
    # Put in cross-stream distance
    latt = lat.copy()
    latt[~np.isfinite(latt)] = -90    
    lont = lon.copy()
    lont[~np.isfinite(lont)] = -90  
    lat1sec = stations[str(sections[secnum][0])].lat
    lat2sec = stations[str(sections[secnum][-1])].lat
    lon1sec = stations[str(sections[secnum][0])].lon
    lon2sec = stations[str(sections[secnum][-1])].lon
    adcpind1 = np.argmin((lat1sec-latt)**2 + (lon1sec - lont)**2)
    adcpind2 = np.argmin((lat2sec-latt)**2 + (lon2sec - lont)**2)
    
    indmax = np.argmax(velmag[adcpind1:adcpind2])
    thetamax = angle[adcpind1+indmax]
    sectionangle = np.arctan2(lat2sec-lat1sec, lon2sec-lon1sec)
    deltaangle = sectionangle - thetamax - np.pi/2
    latmax = latt[adcpind1+indmax]
    lonmax = lont[adcpind1+indmax]
    #normaldistance = distance*np.cos(deltaangle)
    #normaldistane = normaldistance - normaldistance[indmax]
    
    
    for i in sections[secnum]:
        npress = stations[str(i)].press.size
        lats = np.concatenate((lats, np.repeat([stations[str(i)].lat], npress)))   
        lons = np.concatenate((lons, np.repeat([stations[str(i)].lon], npress))) 
        thetas = np.concatenate((thetas, stations[str(i)].theta)) 
        nitrates = np.concatenate((nitrates, stations[str(i)].nitrate))
        phosphates = np.concatenate((phosphates, stations[str(i)].phosphate))
        oxygens = np.concatenate((oxygens, stations[str(i)].oxygen))
    
        rhosort, presssort = stations[str(i)].getRhoSorted()
        press = np.concatenate((press, presssort))   
        rho = np.concatenate((rho, rhosort))
        distance = np.concatenate((distance, np.repeat([stations[str(i)].getDistance(latmax, lonmax)], npress)))*np.cos(deltaangle) # in km
    
    
    return distance, nitrates, phosphates, press, rho
#%%
fig = plt.figure(figsize=(10,6))
axbot = plt.subplot2grid((3,4), (2,0), colspan=4)

latfac = 1e3
    
xi = np.linspace(-100,100,500)*latfac
zi = -np.linspace(-100, 100,200)

#plt.figure()
alldist = []
allpstar = []
for secnum in ['1', '2', '3','4']:
    distance, nitrates, phosphates, press, rho = loadSection(secnum)
    pltvar = (phosphates- nitrates/16)*rho/1000 # mmol per m^3
    mask = (pltvar>-9) & np.isfinite(pltvar)
    maskpress = (press<200) &mask &  (np.abs(pltvar)<0.2)
#    alldist = np.concatenate((alldist, distance[maskpress]))
#    allpstar = np.concatenate((allpstar, pltvar[maskpress]))
    axbot.plot(distance[maskpress], pltvar[maskpress], linestyle='None', marker='x', color='0.8')
    maskpress = (press<800) & (press>200)&mask &  (np.abs(pltvar)<0.2)
    axbot.plot(distance[maskpress], pltvar[maskpress], linestyle='None', marker='+', color='0.8')
    distset = set(distance)
    for loc in distset:
        inds = np.where(distance==loc)
        temppltvar = pltvar[inds]
        temppress = press[inds]
        mask = (temppltvar>-9) & (np.isfinite(temppltvar)) & (temppress<200) & (np.abs(temppltvar)<0.2)
#        plt.plot(loc, np.mean(temppltvar[mask]), marker='x', color='r', linestyle='None')
        alldist = np.concatenate((alldist, [loc]))
        allpstar = np.concatenate((allpstar, [np.mean(temppltvar[mask])]))
        mask = (temppltvar>-9) & (np.isfinite(temppltvar)) & (temppress>200) & (temppress<800) & (np.abs(temppltvar)<0.2)
#        plt.plot(loc, np.mean(temppltvar[mask]), marker='+', color='b', linestyle='None')
        
for secnum in ['1', '2', '3','4']:
    distance, nitrates, phosphates, press, rho = loadSection(secnum)
    pltvar = (phosphates- nitrates/16)*rho/1000 # mmol per m^3
    mask = (pltvar>-9) & np.isfinite(pltvar)
    maskpress = (press<200) &mask &  (np.abs(pltvar)<0.2)
#    alldist = np.concatenate((alldist, distance[maskpress]))
#    allpstar = np.concatenate((allpstar, pltvar[maskpress]))
#    axbot.plot(distance[maskpress], pltvar[maskpress], linestyle='None', marker='x', color='0.8')
    maskpress = (press<800) & (press>200)&mask &  (np.abs(pltvar)<0.2)
#    axbot.plot(distance[maskpress], pltvar[maskpress], linestyle='None', marker='+', color='0.8')
    distset = set(distance)
    for loc in distset:
        inds = np.where(distance==loc)
        temppltvar = pltvar[inds]
        temppress = press[inds]
        mask = (temppltvar>-9) & (np.isfinite(temppltvar)) & (temppress<200) & (np.abs(temppltvar)<0.2)
        plt.plot(loc, np.mean(temppltvar[mask]), marker='x', color='r', linestyle='None')
        alldist = np.concatenate((alldist, [loc]))
        allpstar = np.concatenate((allpstar, [np.mean(temppltvar[mask])]))
        mask = (temppltvar>-9) & (np.isfinite(temppltvar)) & (temppress>200) & (temppress<800) & (np.abs(temppltvar)<0.2)
        plt.plot(loc, np.mean(temppltvar[mask]), marker='+', color='b', linestyle='None')
#    vi = griddata((distance[mask]*latfac, -press[mask]), pltvar[mask], (xi[None,:], zi[:,None]), method='linear', rescale=False)
#    vi[~np.isfinite(vi)] = 0
#    pltvarmld = np.mean(vi[-20:,:],axis=0)
#    pltvardeep = np.mean(vi[-80:-20,:], axis=0)
#    plt.plot(xi, pltvarmld)
#    plt.plot(xi, pltvardeep, linestyle = 'dashed')
plt.ylabel('P* [mmol m$^{-3}$]')
plt.xlabel('Cross-stream Distance')
fitdist = np.linspace(-50, 50, 200)

#fit tanh profile
def curfFunc(x, a, b, amplitude):
    return amplitude*0.5*(1+np.tanh( (x + a)/b))
#def curfFunc(x, a,b,c,d,e):
#    return a*x**4 + b*x**3 + c*x**2 + d*x + e
#def curfFunc(x, a,b,c,d):
#    return a*x**3 + b*x**2 + c*x + d
#popt, pcov = curve_fit(curfFunc, alldist, allpstar, [-12, 10, 0.05])
mask = np.abs(allpstar)<100
popt, pcov = curve_fit(curfFunc, alldist[mask], allpstar[mask], p0 = [20, 2.5, 0.1])

#plt.figure()
#plt.plot(alldist, allpstar, linestyle='None', marker='x')
ct = 13
fitdist = np.linspace(ct - 5, ct + 6, 200)
#plt.plot(fitdist, curfFunc(fitdist, popt[0], 2.5, popt[2]), linestyle='dashed', color='k')
#plt.plot(fitdist, curfFunc(fitdist, popt[0], popt[1], popt[2], popt[3]), linestyle='dashed')
#plt.plot(fitdist, (fitdist-ct)*0.05/5, linestyle='dashed', color='k')

maxflux = 200*(popt[2]*0.5)*1/(2.5*1e3)*100 #d/dx tanh((x+a)/b) = 1/b*sech()^2 max is 1/b


fluxconvpalter = 25*100*0.05/((5e3)**2)*86400*365*1e-3 #kappa * MLD * Delta P / L^2 * sec/year *mol/mmol
plt.ylim((-0.15, 0.25))
plt.grid()
#%% SALINITY PLOT
#axSal = plt.subplot2grid((3,4), (0,2), colspan=2, rowspan=2)
#axSal.set_title('Salinity')
#
#sectdistA, zA, salA, rhoA, yeardayA = ReadAtlantisSection(7)
#
#cmapsal = 'plasma'
#conts = np.linspace(35.25, 36.75, 50)
#im = axSal.contourf(sectdistA, zA, salA, conts, extend='both',cmap=cmapsal)
#for c in im.collections:
#    c.set_edgecolor("face")
#cb = plt.colorbar(im, ax=axSal)
#cb.set_ticks((conts[0], 36, conts[-1]))
#cb.solids.set_edgecolor("face")
#
#axSal.contour(sectdistA, zA, rhoA, 20, colors='k')
#axSal.set_ylim((0, 200))
#axSal.invert_yaxis()
#axSal.set_ylabel('z [m]')
#axSal.set_xlabel('Cross-stream distance [km]')
#%% PHOSPHATE PLOT

cl = 0.1

axPho = plt.subplot2grid((3,4), (0,0), colspan=2, rowspan=2)

distance, nitrates, phosphates, press, rho = loadSection('3')


latfac = 1e3
xi = np.linspace(np.min(distance),np.max(distance),500)*latfac
zi = np.linspace(np.max(press), np.min(press),200)
#zi = -np.linspace(0, 200, 100)
# grid the data.

pltvar = (phosphates- nitrates/16)*rho/1000 # mmol per m^3
#pltvar = oxygens
#pltvar = nitrates
#pltvar = pltvar*rho/np.mean(rho)
mask = (pltvar>-9) & np.isfinite(pltvar)
vi = griddata((distance[mask]*latfac, press[mask]), pltvar[mask], (xi[None,:], zi[:,None]), method='linear', rescale=False)
#vi = griddata(lats[mask], -press[mask],pltvar[mask], xi, zi, interp='linear')
mask = (rho>-9) & np.isfinite(rho)
rhoi = griddata((distance[mask]*latfac, press[mask]), rho[mask], (xi[None,:], zi[:,None]), method='linear', rescale=False)
#rhoi = griddata(lats[mask], -press[mask], rho[mask], xi, zi, interp='linear')
# contour the gridded data, plotting dots at the randomly spaced data points.
cmap = 'plasma'#'plasma'
conts = np.linspace(-cl, cl, 50)
#plt.figure()
im =axPho.contourf(xi/1000, zi, vi, conts, cmap=cmap, extend='both')
for c in im.collections:
    c.set_edgecolor("face")
axPho.set_ylim((0, 1000))
axPho.contour(xi/1000, zi, rhoi, 20, colors='k')
axPho.plot(distance*latfac/1000, press, marker='.', markersize=2,linestyle='None', color='w')

axPho.set_xlabel('Cross-stream distance [km]')
axPho.set_ylabel('Pressure [db]')
axPho.invert_yaxis()
#axPho.set_title('$P^* = PO^{3-}_{4} - \\frac{NO^{-}_{3}}{16}$')
axPho.set_title('Section 3')
plt.subplots_adjust(wspace=1, hspace=0.5)


#%% PHOSPHATE PLOT
axPho = plt.subplot2grid((3,4), (0,2), colspan=2, rowspan=2)

distance, nitrates, phosphates, press, rho = loadSection('4')


latfac = 1e3
xi = np.linspace(np.min(distance),np.max(distance),500)*latfac
zi = np.linspace(np.max(press), np.min(press),200)
#zi = -np.linspace(0, 200, 100)
# grid the data.

pltvar = (phosphates- nitrates/16)*rho/1000 # mmol per m^3
#pltvar = oxygens
#pltvar = nitrates
#pltvar = pltvar*rho/np.mean(rho)
mask = (pltvar>-9) & np.isfinite(pltvar)
vi = griddata((distance[mask]*latfac, press[mask]), pltvar[mask], (xi[None,:], zi[:,None]), method='linear', rescale=False)
#vi = griddata(lats[mask], -press[mask],pltvar[mask], xi, zi, interp='linear')
mask = (rho>-9) & np.isfinite(rho)
rhoi = griddata((distance[mask]*latfac, press[mask]), rho[mask], (xi[None,:], zi[:,None]), method='linear', rescale=False)
#rhoi = griddata(lats[mask], -press[mask], rho[mask], xi, zi, interp='linear')
# contour the gridded data, plotting dots at the randomly spaced data points.
cmap = 'plasma'#'plasma'
#conts = np.linspace(-0.08, 0.08, 50)
#plt.figure()
im =axPho.contourf(xi/1000, zi, vi, conts, cmap=cmap, extend='both')
for c in im.collections:
    c.set_edgecolor("face")
axPho.set_ylim((0, 1000))
axPho.contour(xi/1000, zi, rhoi, 20, colors='k')
axPho.plot(distance*latfac/1000, press, marker='.', markersize=2,linestyle='None', color='w')
#cb = plt.colorbar(im,ax=axPho)
#cb.set_label('mmol m$^{-3}$')
##cb.set_ticks([-0.08, 0, 0.08, 0.16])
#cb.solids.set_edgecolor("face")

axPho.set_xlabel('Cross-stream distance [km]')
axPho.set_ylabel('Pressure [db]')
axPho.invert_yaxis()
# Now adding the colorbar


#axPho.set_title('$P^* = PO^{3-}_{4} - \\frac{NO^{-}_{3}}{16}$')
axPho.set_title('Section 4')
plt.subplots_adjust(wspace=1, hspace=0.5, right=0.85)
cbaxes = fig.add_axes([0.87, 0.41, 0.02, 0.45]) 
cb = plt.colorbar(im, cax = cbaxes)  
cb.set_ticks([-cl, 0, cl])
cb.set_label('$P^* = PO^{3-}_{4} - \\frac{NO^{-}_{3}}{16}$\n[mmol m$^{-3}$]')
#plt.tight_layout()2#%%
#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/FigurePStar.pdf', bbox_inches='tight')
#%%
#vl = 1
#plt.figure()
#plt.pcolor(xi, zi, np.gradient(vi, axis=-1), vmin=-vl, vmax=vl)
#plt.ylim((-750, 0))
#plt.colorbar()
#plt.contour(xi, zi, rhoi, 20, colors='k')
#plt.clim((-1e-3, 1e-3))