#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:47:10 2024

@author: deems
"""

from hebraud_lequeux import hebraud_lequeux
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#Fixed parameters
Gpass = 1.
Gact = 1.
gammacrit = 1.
tau = 1.
cstresspass = Gpass*gammacrit

#Varied parameters
yieldratio = 1.2
birthfrac = 0.625
cstressact = yieldratio*cstresspass
alphacrit = 0.5*birthfrac*cstresspass**2 + 0.5*(1-birthfrac)*cstressact**2
Gcrit = birthfrac/(1. - birthfrac)

alphavec = np.array([0.2,0.4,0.6,0.8])
shearratevec = np.logspace(-5,-1,3000)
stressvec = np.empty(len(shearratevec))

cmap = mpl.colormaps['jet']
colors = cmap(np.linspace(0, 0.5, len(alphavec))) #Don't need to use the full colour map for a small number of values
plt.figure()

for j in range(len(alphavec)):
    
    trialx = 0.5 #inital guess for fsolve

    alpha = alphavec[j] 
    stress_large = (birthfrac*Gpass - (1.-birthfrac)*Gact)*shearratevec #analytic expression for stress at large shearrate
    
    for i in range(0,len(shearratevec)):
        
        shearrate = shearratevec[i]
    
        stressvec[i], trialxcurr = hebraud_lequeux(Gpass, Gact, gammacrit, yieldratio, birthfrac, alpha, tau, shearrate, trialx).meanstress()
        
        trialx = trialxcurr #update the initial guess for fsolve with the previous value's solution
      
    plt.plot(shearratevec,stressvec,linewidth=3,color=colors[j],label=f'alpha = {alpha:.3f}',linestyle='-')
    #plt.plot(shearratevec,stress_large,linewidth=3,color=colors[j],label='large approx',linestyle='-.') 
    #Maybe write a custom log log function that can handle negative values of stress

figname = f'./plots/curves_stress_sherarrate_f={birthfrac}.pdf'
plt.xlabel('shear rate')
plt.ylabel('stress')
plt.title(f'f = {birthfrac}, yieldratio = {yieldratio} \n Gcrit = {Gcrit:.3f}, alpha_crit = {alphacrit:.3f}')
plt.axhline(0,color='black',linestyle='-.') # stress = 0
plt.legend()
plt.savefig(figname)

