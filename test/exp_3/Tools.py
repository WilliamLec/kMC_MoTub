#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 2019
@author: William Lecompte 
Description: 
Library of function for making measurments during the kMC simulation 
of the dynamics of kinesin-1 and MT  
"""


import numpy as np


def measurement_motor(lat):
    Ne = 0 
    Nk = 0
    for i in range(len(lat[0,:,0])):
        for j in range(len(lat[0,0,:])):
            if lat[0,i,j] == 3 or lat[0,i,j] == 4:
                Ne += 1 
            elif lat[1,i,j] == 1:
                Nk += 1
    return (Nk,Ne)
#-----------------------------------------------------------------------------
#--------------Motor flux measurement ---------------
def flux(lat,reac,x,J):
    x = int(x)
    if reac == 2 : # if the next reaction is a walking reaction 
        J[x] = J[x] + 1 # Add 1 to the number of motors which passed the site at position x 
    return J
    
def dens(kp,km,k0):
    rho = ((4*kp+km)-np.sqrt(km**2+4*kp*(km+k0)))/(2*(km+4*kp-k0))
    return rho 
    
def attachment_rate_cte(rho,km,k0):
    kp = (km*rho*(1-rho)-k0*rho**2)/(1-2*rho)**2
    return kp
    
