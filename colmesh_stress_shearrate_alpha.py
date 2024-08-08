
from getdata import getdata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#Folder name where the data is stored is labelled according to the parameters of log spacing 
lowerlim = -5
upperlim = 2
resolution = 1401
#Look at data within a certain range of shear rate
minidx = 0
maxidx = 1201 

#Fixed parameters
Gpass = 1.
Gact = 1.
tau = 1.
cstresspass = 1.

#Varied parameters
yieldratio = 1.2
cstressact = yieldratio*cstresspass
birthfrac = 0.625
alphacrit = 0.5*birthfrac*cstresspass**2 + 0.5*(1-birthfrac)*cstressact**2

#alphavec = np.linspace(0.1,0.9,33)
#alphavec = np.linspace(0.25,1.5,51)
alphavec = np.linspace(0.25,1.,31)
stress = np.empty([len(alphavec),maxidx])

for i in range(len(alphavec)):
    
    alpha = alphavec[i]
    
    stressvec = getdata(Gact, yieldratio, birthfrac, alpha, lowerlim, upperlim, resolution).stress(minidx,maxidx)
    shearratevec = getdata(Gact, yieldratio, birthfrac, alpha, lowerlim, upperlim, resolution).shearrate(minidx,maxidx)
    
    stress[i,:] = stressvec
    
figname = f'./plots/stress_colour_alpha_G2={Gact}_f={birthfrac}_y={yieldratio}.pdf'
plt.figure()
    
col_scale = 1.05 #Limit for the colour bar is +/-col_scale*max(abs(values))
plt.pcolormesh(np.log10(shearratevec),alphavec,stress,norm=colors.SymLogNorm(linthresh=0.01,vmin=-col_scale*np.max(np.abs(stress)),vmax=col_scale*np.max(np.abs(stress))),cmap=plt.cm.seismic,rasterized=True)
#linthresh is the region around zero that is treated linearly
plt.colorbar()
#plt.clim(-col_scale*np.max(np.abs(stress)),col_scale*np.max(np.abs(stress)))
plt.clim()
plt.contour(np.log10(shearratevec),alphavec,stress,levels=[0],linewidths=2.,linestyles='--',cmap='Accent')
plt.axhline(alphacrit,color='black',linestyle='-')
plt.title(f'stress, |G2| = {Gact}, f = {birthfrac}, yieldratio = {yieldratio} \n Gcrit = {birthfrac/(1-birthfrac):.2f}, alpha_c = {alphacrit:.2f}')
plt.xlabel('log(shearrate)')
plt.ylabel('alpha')

plt.savefig(figname)

