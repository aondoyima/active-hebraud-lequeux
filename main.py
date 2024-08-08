
from hebraud_lequeux import hebraud_lequeux
import numpy as np
import argparse
import os
import pickle

#Using argparse means we can easily change parameters when we run the code in comand line,
#e.g python main.py -f 0.6. If no options are specified, the default values are used.
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gact", type = float, default = 1., help = "active shear modulus")
parser.add_argument("-t", "--alpha", type = float, default = 0.7, help = "temperature-like parameter")
parser.add_argument("-f", "--birthfrac", type = float, default = 0.5, help = "ratio of active to passice")
parser.add_argument("-y", "--yieldratio", type = float, default = 1.2, help = "active yield stress/passive yield stress")
args = parser.parse_args()

#Varied parameters
Gact = args.gact
alpha = args.alpha
birthfrac = args.birthfrac
yieldratio = args.yieldratio

#Fixed parameters
Gpass = 1.
tau = 1.
gammacrit = 1.

#Shear rate vector - save each bit of data as a pair of vectors with the stress and shear rate.
#We will use log spacing for the shear rate vector.
lowerlim = -5
upperlim = 2
resolution = 1401
shearratevec = np.logspace(lowerlim,upperlim,resolution)

stressvec = np.empty(len(shearratevec))
trialx = 0.5 #Initial guess for fsolve

#Create a folder to save to save data 
dirname = f'./data_logspaceshearrate_{lowerlim}_{upperlim}_{resolution}/yieldratio={yieldratio:.3f}/alpha={alpha:.3f}/f={birthfrac:.3f}/Gact={Gact:.3f}/'

#Make the folders - exist_ok=True means that if a directory already exists, it won't produce and error. 
#Instead it will make a directory at the first point in the nest where there isn't a directory.
#e.g if I do os.makedirs('A/B/C', exist_ok=True), and A/B exists, then C will be made inside B. If none of them exist, A/B/C will be made, and so on.
os.makedirs(dirname, exist_ok=True)

for i in range(0,len(shearratevec)):
        
    shearrate = shearratevec[i]
    
    stressvec[i], trialxcurr = hebraud_lequeux(Gpass, Gact, gammacrit, yieldratio, birthfrac, alpha, tau, shearrate, trialx).meanstress()
        
    trialx = trialxcurr #update the initial guess for fsolve with the previous value's solution
    
#Derivative of stress with respect to shearrate at very small shear rate.
delta = 5e-9
stress_shearrateminusdelta = hebraud_lequeux(Gpass, Gact, gammacrit, yieldratio, birthfrac, alpha, tau, shearratevec[0] - delta, trialx).meanstress()[0]
stress_shearrateplusdelta = hebraud_lequeux(Gpass, Gact, gammacrit, yieldratio, birthfrac, alpha, tau, shearratevec[0] + delta, trialx).meanstress()[0]
#Use finite difference approximation accurate to delta^2
gradientstress = (stress_shearrateplusdelta - stress_shearrateminusdelta)/(2.*delta)

#Save parameters and data in pickle files.
params = {'passive_modulus': Gpass, 'active_modulus': Gact, 'alpha': alpha, 'birthfrac': birthfrac, 'yieldratio': yieldratio, 'tau': tau, 'passive_yield_strain': gammacrit}

pickle.dump(params, open(dirname+'/params.p', 'wb'))

data = {'stress': stressvec, 'shearrate': shearratevec, 'gradient': gradientstress}

pickle.dump(data, open(dirname+'/data.pickle', 'wb'))


