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
def ReadKnorr_ALL():
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
    for secnum in range(1, 71):
        if secnum<10:
            filename = '/home/jacob/dedalus/LATMIX/KNORR SECTIONS/S4sec00{}.mat'.format(str(secnum)) # Front run
        else:
            filename = '/home/jacob/dedalus/LATMIX/KNORR SECTIONS/S4sec0{}.mat'.format(str(secnum)) # Front run

        if secnum==20:
            print('here')
        else:
            matfile = spio.loadmat(filename,struct_as_record=True, squeeze_me=True)
        
    #    print("Keys: %s" % matfile.keys())
    #    keys = list(matfile.keys())[:]
        
            cgrid = matfile['section']
    
            lat = np.concatenate((lat,cgrid['lat'].tolist()))
            lon = np.concatenate((lon, cgrid['lon'].tolist()))
            tmp = np.concatenate((tmp, cgrid['time'].tolist()))
 

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
