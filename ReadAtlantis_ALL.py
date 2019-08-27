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
def ReadAtlantis_ALL():
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

    #%%
    lat = []
    lon = []
    tmp = []
    for secnum in range(1, 16):
        if secnum<10:
            filename = '/home/jacob/dedalus/LATMIX/surveyD_files/SurveyD0{}.mat'.format(str(secnum)) # Front run
        else:
            filename = '/home/jacob/dedalus/LATMIX/surveyD_files/SurveyD{}.mat'.format(str(secnum)) # Front run

        matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)
    
#    print("Keys: %s" % matfile.keys())
#    keys = list(matfile.keys())[:]
    
        cgrid = matfile['cgrid']
   
        lat = np.concatenate((lat,cgrid.lat))
        lon = np.concatenate((lon, cgrid.lon))
        tmp = np.concatenate((tmp, cgrid.time))
 

    return lat, lon, tmp
#%%
#conts = np.linspace(35, 37, 30)
#plt.figure()
#im = plt.contourf(sectdist, z, sal, conts, cmap=cmo.haline)
#cb = plt.colorbar(im)
#cb.set_ticks((conts[0], 36, conts[-1]))
#plt.contour(sectdist, z, rho, 20, colors='k')
#plt.ylim((0, 200))
#plt.gca().invert_yaxis()
