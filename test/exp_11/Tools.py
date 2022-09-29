#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
    
### -------- Lattice Geometry ------------------

#### without seam structure 
def above(lat,i,j):     
    l = len(lat[:,0])
    return (int((i+1)%l),int(j))

def below(lat,i,j):     
    l = len(lat[:,0])  
    return (int((i-1)%l),int(j))

#### with seam structure as S3_start
def left_s3(lat,i,j):     
    l = len(lat[:,0])  
    if j == 0:  
        I1 = (i+1)%l
        I2 = (i+2)%l
        J = len(lat[0,:]) - 1 
    else:
        I1 = i
        I2 = i
        J = j-1     
    return (int(J),int(I1),int(I2))

def right_s3(lat,i,j):
    l = len(lat[:,0])     
    if j == len(lat[0,:]) - 1: 
        I1 = (i-1)%l
        I2 = (i-2)%l
        J = 0 
    else: 
        I1 = i
        I2 = i 
        J = j+1
    return (int(J),int(I1),int(I2))

def ngbr_s3(i,j,lat): 
    a = 0 
    at = 0 
    b = 0 
    bt = 0 
    U = above(lat,i,j)
    D = below(lat,i,j)
    L = left_s3(lat,i,j)
    L1 = [L[1],L[0]]
    L2 = [L[2],L[0]]
    R = right_s3(lat,i,j)
    R1 = [R[1],R[0]]
    R2 = [R[2],R[0]]
    
    if lat[i,j] > 0: # not a GTP-dimer
        if lat[U[0],U[1]] != 0: 
            a += 1
        if lat[D[0],D[1]] != 0:
            a += 1        
            
        if j == 0: 
            if lat[L1[0],L1[1]] != 0:
                b += 1/2
            if lat[L2[0],L2[1]] != 0:
                b += 1/2
            if lat[R1[0],R1[1]] != 0:
                b += 1
                
        elif j == len(lat[0,:]) - 1: 
            if lat[R1[0],R1[1]] != 0:
                b += 1/2
            if lat[R2[0],R2[1]] != 0:
                b += 1/2
            if lat[L1[0],L1[1]] != 0:
                b += 1
        else: 
            if lat[R1[0],R1[1]] != 0:
                b += 1
            if lat[L1[0],L1[1]] != 0:
                b += 1
    else: 
        if lat[U[0],U[1]] > 0:  
            a += 1
        elif lat[U[0],U[1]] < 0:
            at += 1
            
        if lat[D[0],D[1]] > 0:  
            a += 1
        elif lat[D[0],D[1]] < 0:
            at += 1
        
        if j == 0: 
            if lat[L1[0],L1[1]] > 0:  
                b += 1/2
            elif lat[L1[0],L1[1]] < 0:
                bt += 1/2
        
            if lat[L2[0],L2[1]] > 0:  
                b += 1/2
            elif lat[L2[0],L2[1]] < 0:
                bt += 1/2
            
            if lat[R1[0],R1[1]] > 0:  
                b += 1
            elif lat[R1[0],R1[1]] < 0:
                bt += 1
                
        elif j == len(lat[0,:]) - 1:
            if lat[R1[0],R1[1]] > 0:  
                b += 1/2
            elif lat[R1[0],R1[1]] < 0:
                bt += 1/2
                
            if lat[R2[0],R2[1]] > 0:  
                b += 1/2
            elif lat[R2[0],R2[1]] < 0:
                bt += 1/2
                
            if lat[L1[0],L1[1]] > 0:  
                b += 1
            elif lat[L1[0],L1[1]] < 0:
                bt += 1
                
        else: 
            if lat[L1[0],L1[1]] > 0:  
                b += 1
            elif lat[L1[0],L1[1]] < 0:
                bt += 1
                
            if lat[R1[0],R1[1]] > 0:  
                b += 1
            elif lat[R1[0],R1[1]] < 0:  
                bt += 1
    return (a,at,b,bt)

def MT_length_measure(lat):
    n = len(lat[:,0]) ; m = len(lat[0,:])
    L_pf = np.zeros(m)
    #kymo = np.zeros(n)
    for j in range(m):
        L_pf[j] = n-1
        bool = True
        k = n-1
        while bool: 
            if k < 3:
                bool = False
            elif lat[k,j] == 0:  # Present of a dimer at the site (k,j)  
                L_pf[j] -= 1
                k -= 1
            else: 
                bool = False
    # for i in range(n):
    #     kymo[i] = lat[i,:].mean()
    return (L_pf)