#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script makes the basic surface flux plot.
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

plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 14
plt.rcParams['contour.negative_linestyle'] = 'solid'
#matplotlib.rcParams.update({'font.size': 22})
#%% LOAD DATA

filename = '/home/jacob/dedalus/LATMIX/LES_COARE35_forcing.mat' # Front run
matfile = spio.loadmat(filename,struct_as_record=False, squeeze_me=True)

# List all groups
yearday = matfile['yearday']
tau_cs = matfile['tau_cs']
tau_d = matfile['tau_d']
qnet = matfile['net_heat_flux']
tmag = np.sqrt(tau_cs**2 + tau_d**2)


fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.plot(yearday, tmag, color=color, linewidth=2)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylabel('Wind-stress magnitude (N m$^{-2}$)', color=color)  # we already handled the x-label with ax1
ax1.set_ylim((0, 1.2))
ax1.set_xlim((64.5, 66))
ax1.set_xlabel('yearday')
ax2 = ax1.twinx()
plt.axvline(64.86, linestyle='dashed', color='k')
color = 'tab:red'
ax2.plot(yearday, qnet, color=color, linewidth=2)
ax2.set_ylabel('Net surface heat flux (W m$^{-2}$)', color=color)  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim((0, 1200))
ax1.grid()