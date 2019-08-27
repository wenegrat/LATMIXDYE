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
f = h5py.File(filename, 'r')

# List all groups
print("Keys: %s" % f.keys())
keys = list(f.keys())[:]
