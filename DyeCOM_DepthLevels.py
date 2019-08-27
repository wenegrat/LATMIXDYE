#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
COM Statistics

Created on Mon Mar 26 09:44:57 2018

@author: jacob
"""

import numpy as np
import h5py
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import scipy.stats as stats
from astropy.stats import circstats
#%% LOAD
def loadKappa(filename):
    f = h5py.File(filename, 'r')

    # List all groups
#    print("Keys: %s" % f.keys())
#    keys = list(f.keys())[:]
    
    
    dye1 = f['dye1']

    time = f['t'][:,0]
    x = f['x'][:,0]
    z = f['z'][:,0]
    z = z[1:]
    
    dye1 = dye1[:,1:,:]
    nt, nz, nx = dye1.shape
    
    # COM Calcs    
    # define circular coordinate system for averaging
    theta = x/np.max(x)*2*np.pi
    xi = np.cos(theta)
    zi = np.sin(theta)
    x2i = np.cos(2*theta)
    z2i = np.sin(2*theta)
    
    # Pick which variable to use
    dyevar = dye1
    
    # Pre-allocate
    ZCM = np.zeros((nt,))
    xposcom = np.zeros((nt,))
    xs = np.zeros((nt, nx))
    XVar = np.zeros((nt,nz))
    
    dyehor = np.zeros((nt, nz))
    
    for j in range(0, nz):
        print(j)
        for i in range(0,nt):
            dyehor[i, j] = integrate.trapz(dyevar[i,j,:], x=x) # integrate horizontally at each depth for later vert COM calc
            # Calculate zonal center of mass
            xii = xi*dyevar[i,j,:] #Weighted unwrapped coordinate system
            zii = zi*dyevar[i,j,:]
            tbar = np.arctan2(-np.mean(zii), -np.mean(xii)) + np.pi
#            vbar = np.arctan2(-np.mean(z2i*dyevar[i,j,:]), -np.mean(x2i*dyevar[i,j,:])) + np.pi
            ZCM[i] = np.max(x)*tbar/(2*np.pi) # This is the horizontal (zonal) COM
            
            # Find index of Zonal Center of Mass
            ind = np.where(x>=ZCM[i])[0][0]
            xposcom[i] = x[ind]
            indshift = ind-nx/2 # This is used to center the dye in the domain.
            indshift = indshift.astype('int64')
            
            #Circular shfit everything
            dyeroll = np.roll(dyevar[i,j,:], -indshift, axis=0) # This puts the COM in the center of the domain
#            xroll = np.roll(x, -indshift)
            xs = x - x[int(nx/2)] # Set center coordinate to 0

            
            # Calculate dye variance around COM
            M20 = integrate.trapz(xs**2*dyeroll, x=xs, axis=0) #
            XVar[i,j] = M20/dyehor[i,j] # maybe this sets x=0 at center? Not important for Kappa Calc....
           
            #dyenorm = dyevar[i,j,:]/np.max(dyevar[i,j,:])
            #XVar[i,j] = stats.circvar(np.arctan2(-zi*dyenorm, -xi*dyenorm))*np.max(x)**2/(2*np.pi)**2
    #        XVar[i,j] = circstats.circvar(theta, weights = dyevar[i,j,:]/np.sum(dyevar[i,j,:]))*np.max(x)**2/(2*np.pi)**2
            #test = circstats.circmoment(theta, p=2.0,centered=False, weights = dyevar[i,j,:]/dyehor[i,j])
            #XVar[i,j] = test[0]*np.max(x)**2

    #Calculate the vertical COM at each timestep
    VCMind = np.zeros((nt,))
    Vsigma = np.zeros((nt,))
    ZCOM = np.zeros((nt,))
    for i in range(0, nt):
       VCOM = integrate.trapz(z*dyehor[i,:], x=z)/integrate.trapz(dyehor[i,:], x=z) # restrict z coordinate
       ind = np.argmin(np.abs(z-VCOM))
       if not ind==0:
           VCMind[i] = np.where(z>=VCOM)[0][0]
           ZCOM[i] = z[int(VCMind[i])]
           zs = z - z[int(VCMind[i])]
           Vsigma[i] = integrate.trapz(zs**2*dyehor[i,:], x=zs)/integrate.trapz(dyehor[i,:], x=zs)
       
              
    # Calculate kappa and buoyancy class
    mask = (np.zeros((nt,nz)))
    zu = np.zeros((nt,))
    zd = np.zeros((nt,))
    zdiff = np.zeros((nt,))
    
    indz = np.argmin(np.abs(z+0))
    for i in range(0, nt):
        # 2 sigma method a la Miles
        indu = np.argmin(np.abs(z - (ZCOM[i] + 2*Vsigma[i]**(1/2))))
        if indu>indz: 
            indu = indz
        indd = np.argmin(np.abs(z - (ZCOM[i] - 2*Vsigma[i]**(1/2))))
        # Pick a specific depth range
        #indu = np.argmin(np.abs(z - (0)))
        #indd = np.argmin(np.abs(z - (-150)))
        # Restrict to percent maximum hor average?
        maxd = np.max(dyehor[i,:])
        cutoff = 0.01*maxd
        
        zu[i] = z[indu]
        zd[i] = z[indd]
        zdiff[i] = z[indu] - z[indd]
        mask[i,indd:indu] = 1
    #    mask[i,(dyehor[i,:]/5e3 < 0.005) ] = 0
    #XVarB = integrate.trapz(XVar*mask, x=z, axis=-1)/np.abs(z[0])
    temp = XVar**(1/2)*mask*dyehor
    temp[np.isnan(temp)] = 0
    #XVarB = (integrate.trapz(temp, x=z, axis=-1)/(integrate.trapz(dyehor*mask, x=z, axis=-1)))**2
    XVarB = integrate.trapz(XVar*mask*dyehor, x=z, axis=-1)/(integrate.trapz(dyehor*mask, x=z, axis=-1))
    
    #XVarB = integrate.trapz(XVar*mask, x=z, axis=-1)/(np.abs(zdiff))
    
    kappa = 0.5*np.gradient(XVarB[0:time.size], axis=0)/np.gradient(time)
    return kappa, XVarB, time

#%% 
    kappa_front, Var_front, time_front = loadKappa('/home/jacob/dedalus/LATMIX/run10k.mat') # Pick which run
    kappa_nofront, Var_nofront, time_nofront = loadKappa('/home/jacob/dedalus/LATMIX/nofront.mat') # Pick which run
    kappa_nogeo, Var_nogeo, time_nogeo = loadKappa('/home/jacob/dedalus/LATMIX/notw.mat') # Pick which run
#%%
indl = np.argmin(np.abs(time_front/86400+64.5 - 64.88))
indr = np.argmin(np.abs(time_front/86400+64.5 - 65.38))

pstargrad = 0.06/2.5e3
#pstargrad = 0.1/1e3
bld = 100
kappaz = kappa_front
#kappaz[kappaz<0] = 0
avgkappa = np.nanmean(kappaz[indl:indr])
peakkappa = np.max(kappaz[indl:indr])
kappasteady = avgkappa*5
kappaoscill = avgkappa/5
maxkappa = np.nanmax(kappaz)

print('K average: ' + str(avgkappa))
print('K peak:    ' + str(peakkappa))

print('K steady: ' + str(kappasteady))
print('K oscill: ' + str(kappaoscill))

#print('Averaged: ' + str(avgkappa*bld*pstargrad)) # estimating one event per 3 days
#print('Maximum:  ' + str(maxkappa*bld*pstargrad)) # estimating one event per 3 days

#%% PLOT KAPPA
plt.figure()
plt.plot(time_front/86400+64.5,kappa_front[0:time_front.size])
plt.plot(time_nofront/86400+64.5,kappa_nofront[0:time_nofront.size])
plt.plot(time_nogeo/86400+64.5,kappa_nogeo[0:time_nogeo.size])

plt.xlabel('Yearday')
plt.ylabel('$\kappa_h$ m$^2$s$^{-1}$')
plt.grid()
plt.tight_layout()
plt.xlim((64.5, 66))

#%% PLOT Variance
plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 20
plt.rcParams['contour.negative_linestyle'] = 'solid'

plt.figure(figsize=(12,4))
plt.plot(time_front/86400+64.5,Var_front[0:time_front.size], label='FRONT', linewidth=2)
plt.plot(time_nofront/86400+64.5,Var_nofront[0:time_nofront.size], label='NO-FRONT', linewidth=2)
plt.plot(time_nogeo/86400+64.5,Var_nogeo[0:time_nogeo.size], label='NO-GEOMIX', linewidth=2)

plt.xlabel('Yearday')
plt.ylabel('Cross-frontal variance, m$^2$')
plt.grid()
plt.tight_layout()
plt.legend()
plt.xlim((64.87, 66))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#plt.savefig('/home/jacob/Dropbox/GulfStreamDye/LATMIXSCIENCE/FigureModelVariance.pdf', bbox_inches='tight')

#%% DO SOME VALIDATION

xg = np.linspace(-20e3, 20e3, 1000)
#time = np.linspace(0, 10, 500)
time = time_front
dye = np.zeros((time.size, xg.size))

sigma2 = np.zeros((time.size,))
sigma2[0] = 75**2
dye[0,:] = np.exp( -(xg**2/(2*25**2) ))

kap = 90
for i in range(1, time.size):
    kap = kappa_front[i]
    if ~np.isfinite(kap):
        kap = 0
    sigma2[i] = sigma2[i-1] + (time[i] - time[i-1])*2*kap
    dye[i,:] = np.exp( -(xg**2/(2*sigma2[i]) ))
plt.figure()
#plt.plot(x, dye[0,:])
tl = range(0, np.argmin(np.abs(time_front/86400-65.5+64.5)), 50)
plt.plot(xg, np.transpose(dye[tl,:]))
plt.legend()

plt.figure()
plt.plot(time_front/86400, sigma2)
plt.plot(time_front/86400, Var_front[0:time_front.size])



#%% PICK A PARTICULAR TIME AND PLOT THE INTEGRATED DYE VARIANCE
# 1) Load a run
# 2) Pick a time
# 3) Calculate vertically integrated dye concentration
# 4) Compare it to the Gaussian shape inferred from above integration

f = h5py.File('/home/jacob/dedalus/LATMIX/run10k.mat', 'r')
dye1 = f['dye1']
time = f['t'][:,0]
x = f['x'][:,0]
z = f['z'][:,0]
z = z[1:]
dye1 = dye1[:,1:,:]
nt, nz, nx = dye1.shape

# COM Calcs    
# define circular coordinate system for averaging
theta = x/np.max(x)*2*np.pi
xi = np.cos(theta)
zi = np.sin(theta)
x2i = np.cos(2*theta)
z2i = np.sin(2*theta)

#%% CONTINUED...
#Pick a time
timetarget = 65.3

timeind = np.argmin(np.abs(time/86400+64.5-timetarget))

dyevert = integrate.trapz(dye1[timeind,:,:], x=z, axis=0)

# Calculate zonal center of mass
xii = xi*dyevert #Weighted unwrapped coordinate system
zii = zi*dyevert
tbar = np.arctan2(-np.mean(zii), -np.mean(xii)) + np.pi
#            vbar = np.arctan2(-np.mean(z2i*dyevar[i,j,:]), -np.mean(x2i*dyevar[i,j,:])) + np.pi
ZCM = np.max(x)*tbar/(2*np.pi) # This is the horizontal (zonal) COM
            
 # Find index of Zonal Center of Mass
ind = np.where(x>=ZCM)[0][0]
xposcom = x[ind]
indshift = ind-nx/2 # This is used to center the dye in the domain.
indshift = indshift.astype('int64')
            
            #Circular shfit everything
dyeroll = np.roll(dyevert, -indshift, axis=0) # This puts the COM in the center of the domain
#            xroll = np.roll(x, -indshift)
xs = x - x[int(nx/2)] # Set center coordinate to 0


vsig = integrate.trapz(xs**2*dyeroll, x=xs)/integrate.trapz(dyeroll, x=xs)
plt.figure()
plt.plot(xs, dyeroll)
plt.plot(xg, dye[timeind,:]*np.max(dyeroll)/np.max(dye[timeind,:]))
plt.axvline(np.sqrt(vsig), color='r')
plt.axvline(np.sqrt(sigma2[timeind]), color='r', linestyle='dashed')

    #%% PLOT VARIANCE
plt.figure()
plt.plot(time/86400+64.5,4*np.sqrt(XVarB[0:time.size]))
plt.xlabel('Yearday')
plt.ylabel('4$\sigma$ [m]')
plt.grid()
plt.tight_layout()
plt.xlim((64.5, 66))
#%%

plt.figure()
plt.pcolor(time/86400,z, np.transpose(kappa), vmin=0, vmax=500)
plt.colorbar()
plt.plot(time/86400, ZCOM[0:1569])

#%%
plt.figure()
plt.pcolor(time/86400 +64.5,z, np.transpose(XVar[0:1569,:]*mask[0:1569,:]))
plt.colorbar()
plt.plot(time/86400 + 64.5, ZCOM[0:1569])
plt.plot(time/86400 + 64.5, zd[0:1569])
plt.plot(time/86400 + 64.5, zu[0:1569])
plt.xlabel('Yearday')
plt.xlim((64.5, 66))
plt.ylabel('z')
#%%
plt.figure()
plt.pcolor(time/86400+64.5,z, np.transpose(dyehor[0:1569,:]))