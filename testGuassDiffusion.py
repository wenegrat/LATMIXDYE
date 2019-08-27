#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:26:18 2019

@author: jacob
"""
import numpy as np
# Test wrapped vs unwrapped

#%%
k = [0, 1, 2, 3, 4, 5]
theta = np.linspace(0, 2*np.pi, 100)
gaussW = np.zeros((theta.shape))
sigma = 20
for i in k:
    gaussW += 1/(np.sqrt(2*np.pi*sigma**2))*np.exp(- (theta - np.pi)**2/(2*sigma**2))
    
    
x = np.linspace(0, 2*np.pi, theta.size)

#theta = x/2*2*np.pi
#xi = np.cos(theta)
#zi = np.sin(theta)


plt.figure()
#plt.plot(x, guass1)
plt.plot(theta, gaussW)


xbar = np.average(x, weights=gaussW)
xvar = np.average((x-xbar)**2, weights=gaussW)

print('Actual sigma^2: ' + str(sigma**2))
print('Inferred sigma^2: ' + str(xvar))