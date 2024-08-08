#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 21:50:01 2024

@author: deems
"""

from getdata import getdata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

lowerlim = -5
upperlim = 2
resolution = 1401

#Fixed parameters
Gpass = 1.
Gact = 1.
tau = 1.
cstresspass = 1.

#Varied parameters
yieldratio = 1.2
cstressact = yieldratio*cstresspass

alphavec = np.linspace(0.25,1.,31)
#alphavec = np.linspace(0.25,1.5,51)
birthfracvec = np.linspace(0.,1.,41)
alphacrit = 0.5*birthfracvec*cstresspass**2 + 0.5*(1-birthfracvec)*cstressact**2

#stress and gradient of the stress at small strain
stress_small = np.empty([len(alphavec),len(birthfracvec)])
gradient_small = np.empty([len(alphavec),len(birthfracvec)]) 

for i in range(len(alphavec)):
    for j in range(len(birthfracvec)):
        
        alpha = alphavec[i]
        birthfrac = birthfracvec[j]
        
        stress_small[i,j] = getdata(Gact, yieldratio, birthfrac, alpha, lowerlim, upperlim, resolution).stress_small() 
        gradient_small[i,j] = getdata(Gact, yieldratio, birthfrac, alpha, lowerlim, upperlim, resolution).gradient()
        #I only calculated the gradiet for the smallest value of shear rate used.


col_scale = 1.05 #Limit for the colour bar is +/-col_scale*max(abs(values))

plt.figure(1)         
colmesh1 = plt.pcolormesh(birthfracvec,alphavec,stress_small,norm=colors.SymLogNorm(linthresh=0.2,vmin=-col_scale*np.max(np.abs(stress_small)),vmax=col_scale*np.max(np.abs(stress_small))),cmap=plt.cm.seismic,rasterized=True)
#The data suggests log spacing is better for the colours 
plt.colorbar(colmesh1)
plt.clim(-col_scale*np.max(np.abs(stress_small)),col_scale*np.max(np.abs(stress_small)))
plt.contour(birthfracvec,alphavec,stress_small,levels=[0],linewidths=2.,linestyles='-',cmap='gist_gray')
plt.axvline(0.5,color='black',linestyle=':') #Marks f = 0.5. We need f > 0.5 for the system to be passive at large shear rate.
plt.title(f'stress at small shear rate, yieldratio = {yieldratio}')
plt.xlabel('f')
plt.ylabel('alpha')
plt.plot(birthfracvec,alphacrit,linewidth=2,color='black',linestyle = '--')

figname1 = f'./plots/stress_low_shearrate_{yieldratio}.pdf'
plt.savefig(figname1)

plt.figure(2)
colmesh2 = plt.pcolormesh(birthfracvec,alphavec,gradient_small,norm=colors.SymLogNorm(linthresh=10.,vmin=-col_scale*np.max(np.abs(gradient_small)),vmax=col_scale*np.max(np.abs(gradient_small))),cmap=plt.cm.seismic,rasterized=True)
#The data suggests log spacing is better for the colours
plt.colorbar(colmesh2)
plt.clim(-col_scale*np.max(np.abs(gradient_small)),col_scale*np.max(np.abs(gradient_small)))
cont2 = plt.contour(birthfracvec,alphavec,gradient_small,levels=[0],linewidths=2.,linestyles='-',cmap='Accent')
cont3 = plt.contour(birthfracvec,alphavec,stress_small,levels=[0],linewidths=2.,linestyles='-',cmap='gist_gray')
plt.axvline(0.5,color='black',linestyle=':') #Marks f = 0.5. We need f > 0.5 for the system to be passive at large shear rate.
plt.title(f'gradient at small shear rate, yieldratio = {yieldratio}')
plt.xlabel('f')
plt.ylabel('alpha')
plt.plot(birthfracvec,alphacrit,linewidth=2,color='black',linestyle = '--')

figname2 = f'./plots/stress_gradient_low_shearrate_{yieldratio}.pdf'
plt.savefig(figname2)
