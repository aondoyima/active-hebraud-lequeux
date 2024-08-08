#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 8 14:40:16 2024

@author: deems
"""

import pickle

class getdata:
    def __init__(self, Gact, yieldratio, birthfrac, alpha, lowerlim, upperlim, resolution):
        
        self.init_params(Gact, yieldratio, birthfrac, alpha, lowerlim, upperlim, resolution)
    
    def init_params(self, Gact, yieldratio, birthfrac, alpha, lowerlim, upperlim, resolution):
        
        dirname = f'./data_logspaceshearrate_{lowerlim}_{upperlim}_{resolution}/yieldratio={yieldratio:.3f}/alpha={alpha:.3f}/f={birthfrac:.3f}/Gact={Gact:.3f}/'
        self.data = pickle.load(open(dirname+'data.pickle', 'rb')) 
        self.params = pickle.load(open(dirname+'params.p', 'rb')) 
        
    def stress(self,minidx,maxidx):
        
        stress = self.data['stress'][minidx:maxidx]
        
        return stress
    
    def shearrate(self,minidx,maxidx):
        
        shearrate = self.data['shearrate'][minidx:maxidx]
        
        return shearrate
    
    def stress_small(self):
        
        stress_small = self.data['stress'][0]
        
        return stress_small
    
    def gradient(self):
        
        #I only calculated the gradient for the smallest value of shear rate used.
        gradient = self.data['gradient']
        
        return gradient
    
    def params(self):
        
        return self.params
    
    
        
        

        
        