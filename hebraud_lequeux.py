
import numpy as np
from scipy.optimize import fsolve

#Here we calculate the 

class hebraud_lequeux:
    def __init__(self, Gpass, Gact, gammacrit, yieldratio, birthfrac, alpha, tau, shearrate, trialx):
        
        self.init_params(Gpass, Gact, gammacrit, yieldratio, birthfrac, alpha, tau, shearrate, trialx)
    
    def init_params(self, Gpass, Gact, gammacrit, yieldratio, birthfrac, alpha, tau, shearrate, trialx):
        
        self.Gpass = Gpass
        self.Gact = Gact
        
        self.gpass = Gpass*shearrate
        self.gact = -Gact*shearrate
        
        self.cstresspass = Gpass*gammacrit
        self.cstressact = yieldratio*self.cstresspass 
        #The above line is equivalent to cstressact = Gact*gammacritact, where gammacritact = yieldratio*Gpass*gammacrit/Gact
                
        self.birthfrac = birthfrac
        self.alpha = alpha
        self.tau = tau
        self.shearrate = shearrate 
        
        self.trialx = trialx
        
    def alphafunc(self,x):
        
        #This is the function that determines the diffusion coefficient via 
        #alphafunc(x) = 0, where x = sqrt(D*tau)
        
        ypass = self.gpass*self.tau/x**2
        xipass = 1. + 4./((x**2)*(ypass**2))
        if self.gpass > 0:
            sqrtxipass = np.sqrt(xipass)
        elif self.gpass < 0:
            sqrtxipass = -np.sqrt(xipass)
        
        yact = self.gact*self.tau/x**2
        xiact = 1. + 4./((x**2)*(yact**2))
        if self.gact > 0:
            sqrtxiact = np.sqrt(xiact)
        elif self.gact < 0:
            sqrtxiact = -np.sqrt(xiact)
        
        fpass = (self.birthfrac*self.cstresspass/ypass)*(1 + (2./(ypass*self.cstresspass) + sqrtxipass)*np.tanh(0.5*ypass*self.cstresspass))/(np.tanh(0.5*ypass*self.cstresspass) + sqrtxipass)
        
        fact = ((1. - self.birthfrac)*self.cstressact/yact)*(1. + (2./(yact*self.cstressact) + sqrtxiact)*np.tanh(0.5*yact*self.cstressact))/(np.tanh(0.5*yact*self.cstressact) + sqrtxiact)
        
        func = x**2 + fpass + fact - self.alpha
        
        return func 
            
    def meanstress(self):
        
        #Solve the nonlinear equation alphafunc(x) = 0, where x = sqrt(D*tau)
        x = fsolve(self.alphafunc, self.trialx)
        D = x**2/self.tau
        
        Gamma = D/self.alpha
        
        ypass = self.gpass/D
        xipass = 1. + 4./(x**2*ypass**2)
        
        #The next 12 lines come from the fact that the probability distributions are exponential and have to decay appropriately 
        #as the magnitude of the shearrate apporaces infinity. Because of this, the sign in front of sqrt(xipass) or sqrt(xiact) 
        #(see supplementary material) is determined by the sign of gpass and gact. 
        if self.gpass > 0:
            sqrtxipass = np.sqrt(xipass)
        elif self.gpass < 0:
            sqrtxipass = -np.sqrt(xipass)
        
        yact = self.gact/D
        xiact = 1. + 4./(x**2*yact**2)
        
        if self.gact > 0:
            sqrtxiact = np.sqrt(xiact)
        elif self.gact < 0:
            sqrtxiact = -np.sqrt(xiact)
        
        apass = x**4*ypass**2 + 2.*x**2 - 2./ypass**2 + 0.5*self.cstresspass**2 + sqrtxipass*self.cstresspass*(x**2*ypass - 1./ypass)
        bpass = (x**4*ypass**2 + 0.5*self.cstresspass**2)*sqrtxipass + x**2*ypass*self.cstresspass + self.cstresspass/ypass
        
        aact = x**4*yact**2 + 2.*x**2 - 2./yact**2 + 0.5*self.cstressact**2 + sqrtxiact*self.cstressact*(x**2*yact - 1./yact)
        bact = (x**4*yact**2 + 0.5*self.cstressact**2)*sqrtxiact + x**2*yact*self.cstressact + self.cstressact/yact
        
        meanstresspass = (self.birthfrac*Gamma/self.gpass)*(np.tanh(0.5*ypass*self.cstresspass)*apass + bpass)/(np.tanh(0.5*ypass*self.cstresspass) + sqrtxipass)
        meanstressact = ((1 - self.birthfrac)*Gamma/self.gact)*(np.tanh(0.5*yact*self.cstressact)*aact + bact)/(np.tanh(0.5*yact*self.cstressact) + sqrtxiact)
        
        meanstress = meanstresspass + meanstressact
        
        return meanstress, x
    
    def probdist(self,s):
        
        x = fsolve(self.alphafunc, self.trialx)
        D = x**2/self.tau
        
        Gamma = D/self.alpha
        
        ypass = self.gpass/D
        xipass = 1. + 4./(x**2*ypass**2)
        if self.gpass > 0:
            sqrtxipass = np.sqrt(xipass)
        elif self.gpass < 0:
            sqrtxipass = -np.sqrt(xipass)
        
        yact = self.gact/D
        xiact = 1. + 4./(x**2*yact**2)
        if self.gact > 0:
            sqrtxiact = np.sqrt(xiact)
        elif self.gact < 0:
            sqrtxiact = -np.sqrt(xiact)
            
        Gammapass = self.birthfrac*Gamma
        Gammaact = (1. - self.birthfrac)*Gamma
        
        Ppass = 0.
        Pact = 0.
        
        Apass = 1./((1. + sqrtxipass)*np.exp(0.5*ypass*self.cstresspass) - (1. - sqrtxipass)*np.exp(-0.5*ypass*self.cstresspass))
        
        if s <= -self.cstresspass:
            Ppass = 2*(Gammapass/self.gpass)*np.exp(0.5*ypass*sqrtxipass*self.cstresspass)*np.exp(0.5*ypass*(1. + sqrtxipass)*s)*Apass
        elif s > -self.cstresspass and s <= 0.:
            Ppass = (Gammapass/self.gpass)*(1. + sqrtxipass)*np.exp(-0.5*ypass*self.cstresspass)*(np.exp(ypass*(s + self.cstresspass)) + (1. - sqrtxipass)/(1. + sqrtxipass))*Apass
        elif s > 0. and s <= self.cstresspass:
            Ppass = (Gammapass/self.gpass)*(1. - sqrtxipass)*np.exp(0.5*ypass*self.cstresspass)*(np.exp(ypass*(s - self.cstresspass)) + (1. + sqrtxipass)/(1. - sqrtxipass))*Apass
        elif s > self.cstresspass:
            Ppass = 2*(Gammapass/self.gpass)*np.exp(0.5*ypass*sqrtxipass*self.cstresspass)*np.exp(0.5*ypass*(1. - sqrtxipass)*s)*Apass
            
        Aact = 1./((1. + sqrtxiact)*np.exp(0.5*yact*self.cstressact) - (1. - sqrtxiact)*np.exp(-0.5*yact*self.cstressact))
        
        if s <= -self.cstressact:
            Pact = 2*(Gammaact/self.gact)*np.exp(0.5*yact*sqrtxiact*self.cstressact)*np.exp(0.5*yact*(1. + sqrtxiact)*s)*Aact
        elif s > -self.cstressact and s <= 0.:
            Pact = (Gammaact/self.gact)*(1. + sqrtxiact)*np.exp(-0.5*yact*self.cstressact)*(np.exp(yact*(s + self.cstressact)) + (1. - sqrtxiact)/(1. + sqrtxiact))*Aact
        elif s > 0. and s <= self.cstressact:
            Pact = (Gammaact/self.gact)*(1. - sqrtxiact)*np.exp(0.5*yact*self.cstressact)*(np.exp(yact*(s - self.cstressact)) + (1. + sqrtxiact)/(1. - sqrtxiact))*Aact
        elif s > self.cstressact:
            Pact = 2*(Gammaact/self.gact)*np.exp(0.5*yact*sqrtxiact*self.cstressact)*np.exp(0.5*yact*(1. - sqrtxiact)*s)*Aact
            
        return Ppass, Pact
        
        