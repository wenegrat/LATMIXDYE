#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 12:46:28 2019

@author: jacob
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:28:10 2019

@author: jacob
"""
def ReadAtlantisSection(secnum):
    import h5py
    import scipy.io as spio
    import matplotlib.pyplot as plt
    import numpy as np
    from cmocean import cm as cmo
    from scipy.interpolate import griddata
    from scipy.optimize import curve_fit
    #from matplotlib.mlab import griddata
    import datetime as dt
    import gsw
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.size'] = 14
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    #%%
    def getDistance(lat, lon, startlat, startlon):
          latdeg = lat*np.pi/180
          londeg = lon*np.pi/180
          startlatdeg = startlat*np.pi/180
          startlondeg = startlon*np.pi/180
          
          d=2*np.arcsin(np.sqrt((np.sin((startlatdeg-latdeg)/2))**2 + 
                     np.cos(startlatdeg)*np.cos(latdeg)*(np.sin((startlondeg-londeg)/2))**2))
          dkm = 2*np.pi*6356.752*d/(2*np.pi)*np.sign(latdeg-startlatdeg)
          return dkm
    def getRho(sal, theta, press, lon, lat):
          SA = gsw.SA_from_SP(sal,press, lon, lat)
          CT = gsw.CT_from_t(SA, theta, press)
          rho = gsw.rho(SA, CT, press)
          return rho
    #%%
    filename = '/home/jacob/dedalus/LATMIX/surveyD_files/SurveyD0{}.mat'.format(str(secnum)) # Front run
    matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)
    
#    print("Keys: %s" % matfile.keys())
#    keys = list(matfile.keys())[:]
    
    cgrid = matfile['cgrid']
    u = cgrid.u
    v = cgrid.v
    sal = cgrid.sal
    theta = cgrid.t
    lat = cgrid.lat
    z = cgrid.depths
    lon = cgrid.lon
    tmp = cgrid.time
    nz, nx = u.shape
    
    rho = np.nan*u
    for i in range(0, nx):
        rho[:,i] = getRho(sal[:,i], theta[:,i], z, lon[i], lat[i])
        
    
    day = dt.datetime.fromordinal(tmp.astype('int')[0])
    dayfrac = dt.timedelta(days=tmp[0]%1) - dt.timedelta(days = 366)
    yearday = day + dayfrac 
    yearday = yearday - dt.datetime(yearday.year, 1, 1, 0, 0)
    
    #print(yearday)
#    yearday = 30+29 + tmp/86400 # NEED TO CHECK THIS, DOESN"T MATCH LEIF"S MATLAB CODE
    
    #%% Calculate direction of current
    vtemp = v.copy()
    vtemp[~np.isfinite(vtemp)] = 0
    utemp = u.copy()
    utemp[~np.isfinite(utemp)] = 0
    ubar = np.mean(utemp, axis=0)
    vbar = np.mean(vtemp, axis=0)
    angle = np.arctan2(vbar, ubar)
    velmag = np.sqrt(ubar**2 + vbar**2)
    
    indmax = np.argmax(velmag)
    thetamax = angle[indmax]
    lat2sec = np.max(lat)
    lat1sec = np.min(lat)
    lon2sec = np.max(lon)
    lon1sec = np.min(lon)
    sectionangle = np.arctan2(lat2sec-lat1sec, lon2sec-lon1sec)
    deltaangle = sectionangle - thetamax 
    latmax = lat[indmax]
    lonmax = lon[indmax]
    
    sectdist = getDistance(lat, lon, latmax, lonmax)*np.cos(deltaangle)

    return sectdist, z, sal, rho, yearday
#%%
#conts = np.linspace(35, 37, 30)
#plt.figure()
#im = plt.contourf(sectdist, z, sal, conts, cmap=cmo.haline)
#cb = plt.colorbar(im)
#cb.set_ticks((conts[0], 36, conts[-1]))
#plt.contour(sectdist, z, rho, 20, colors='k')
#plt.ylim((0, 200))
#plt.gca().invert_yaxis()
