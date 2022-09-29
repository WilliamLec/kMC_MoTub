"""
Created on May 202
@author: William Lecompte
Description: 
kMC simulation based on the VSS method for modelling 
a simple 2-head motors walk along a 2D lattice + local scan + list of list
"""
import numpy as np
import random as rd
import time
import kMC_Library as kMC
import Tools as tl

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("seed", type=float, help="seed")
parser.add_argument("L", type=int, help="length")
parser.add_argument("ite", type=float, help="ite")
args = parser.parse_args()

ite = args.ite
rd.seed(args.seed) # Seed initialization 
tick_max = 10**4 # Stop condiontion on loops 
N_stat = 10 
data = np.zeros((N_stat,2))

""" Parameters for initial lattice configuration 
     and for the motor walk model 
"""

n = args.L # Nbr of dimer 
m = 1 # Nbr of Protofilaments 

Nr = 3 # Nbr of rate constants
kw = 100 # Cte rate of walking process  
km = 1 # Cte rate of detachment process 
rho = 0.01 # motor density
kp = tl.attachment_rate_cte(rho, km, 0) # Cte rate of attachment process 

""" loop over N_stat time for statistics purposes
"""
for i in range(N_stat):    
    """ Initialize the lattice + event lists 
        lat: lattice of size (n,m)
        k: rate constant vector
        occ: number of enabled site for each reaction type, it's the occurence vector
        ersl[j] = list of the enable site positions for reaction j 
        lsre[k] = inverse of ersl = list of the reaction enable for site with index k 
    """
    t_init = time.time()
    t = 0 
    tick = 0
    lat = kMC.latice_init(n, m)
    k = kMC.rate_constant_init(Nr, kp, kw, km)
    occ, ersl, lsre = kMC.initialize_eventlist(lat, Nr)
    
    
    """ loop the simulation while the stop condition is not verified
    """
    while tick < tick_max: 
        r = rd.random() # 1st draw 
        C, reac = kMC.next_reaction_draw(k, occ, r)
        
        r = rd.random() # 2nd draw
        t += -(np.log(r)/C) # Simulated time evolution
        
        r = rd.randint(0,int(occ[reac])-1) # 3rd draw
        loc = ersl[reac][r]
        
        lat, occ, ersl, lsre = kMC.update_lattice_configuration(lat, occ, ersl, lsre, reac, loc)
        tick += 1
        """ Monitor the simulation progression every 500 loops
        """
        if tick%500 == 0: 
            compu_t = time.time()-t_init
            out = open('lattice_'+'%.f'%n+'_set_'+'%.f'%ite+'_output.txt','w')
            out.write('loop='+'%.f'%tick+', ite='+'%.2f'%((i/N_stat)*100)+', t_simu ='+'%.2f'%(compu_t/60)+' mins \n')
            out.close()
        
    
    """ measure the computation time to compute the loop velocity
    """
    t_comput = time.time() - t_init # time to compute the script
    data[i,0] = t_comput/tick # velocity loop 
    data[i,1] = t/tick # mean time step 

np.savetxt('lattice_'+'%.f'%n+'_set_'+'%.f'%ite+"_efficiency_measurement", data)