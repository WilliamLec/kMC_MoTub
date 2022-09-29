"""
Created on April 2019
@author: William Lecompte
Description: 
kMC simulation based on the VSS method for modelling 
a simple 2-head motors walk along a 2D lattice
"""
import numpy as np
import random as rd
import time
import matplotlib.pyplot as plt
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
data = np.zeros((N_stat, 2))

""" Parameters for initial lattice configuration 
     and for the MoTub model 
"""
BOX_LENGTH = args.L
BOX_WIDTH = 10
LATTICE_HEIGHT = int(0.5*BOX_LENGTH)
SEED_HEIGHT = 3
CAP_HEIGHT = 3

DEFECT_X_POSITION = int(0.5*LATTICE_HEIGHT)
DEFECT_Y_POSITION = int(0.5*BOX_WIDTH)

Nr = 7 + (3*3*5*5)*2 # Nbr of rate constants
k_MOTOR_DETACH = 1 # in s-1 
ENHANCED_FACTOR = 100  
MOTOR_DENSITY = 0.01 # motor density
k_MOTOR_ATTACH = tl.attachment_rate_cte(MOTOR_DENSITY, k_MOTOR_DETACH, 0) 
k_WALK = 100 # in s-1

HYDROLYSIS_TIME = 5/3 # in s
k_DIMER_ATTACH = 20 # in s-1
RELAXATION_TIME = 0.1 # in s

TOTAL_BINDING_ENERGY_GDP_TUBULIN = 45 # in kT
ANISOTROPY = 0.5
LONGITUDINAL_ENERGY_GDP_TUBULIN = (1/(2*(1+ANISOTROPY)))*TOTAL_BINDING_ENERGY_GDP_TUBULIN
LATERAL_ENERGY_GDP_TUBULIN = ANISOTROPY*LONGITUDINAL_ENERGY_GDP_TUBULIN

CONSTRAST = 1.28 
STABILIZATION_GTP_TUBULIN = 0.5*(1-CONSTRAST)*TOTAL_BINDING_ENERGY_GDP_TUBULIN
LONGITUDINAL_ENERGY_GTP_TUBULIN = LONGITUDINAL_ENERGY_GDP_TUBULIN + STABILIZATION_GTP_TUBULIN
LATERAL_ENERGY_GTP_TUBULIN = LATERAL_ENERGY_GDP_TUBULIN

MOTOR_DESTABILIZATION = 2 # in kT
DIRECTIONALITY = 1 

""" loop over N_stat time for statistics purposes
"""
for i in range(N_stat):
    """ Initialize the lattice + event lists 
        lat: lattice of size (n,m)
        k: rate constant vector
        occ: number of enabled site for each reaction type, it's the occurence vector
        ersl: enabled reaction site list, it's an one-hot encoding matrice     
    """
    t = 0 # Initial time
    t_init = time.time()
    lat = kMC.latice_initialization(BOX_LENGTH, BOX_WIDTH, LATTICE_HEIGHT, SEED_HEIGHT, CAP_HEIGHT,
                                    DEFECT_X_POSITION, DEFECT_Y_POSITION)
    k = kMC.rate_constant_initizalition(k_MOTOR_ATTACH, k_MOTOR_DETACH, ENHANCED_FACTOR, k_WALK,
                                        HYDROLYSIS_TIME, k_DIMER_ATTACH, RELAXATION_TIME, LONGITUDINAL_ENERGY_GDP_TUBULIN,
                                        LATERAL_ENERGY_GDP_TUBULIN, LONGITUDINAL_ENERGY_GTP_TUBULIN,
                                        LATERAL_ENERGY_GTP_TUBULIN, MOTOR_DESTABILIZATION)
    occ, ersl = kMC.initialize_event_lists(lat, Nr, DIRECTIONALITY)
    
    """ Initialized observables 
    """
    tick = 0 # number of loops 

    
    """ loop the simulation while the stop condition is not verified
    """
    while tick < tick_max: 
        r = rd.random() # 1st draw 
        C, reac = kMC.next_reaction_draw(k, occ, r)
        
        r = rd.random() # 2nd draw
        t += -(np.log(r)/C) # Simulated time evolution
        
        r = rd.randint(0,int(occ[reac])-1) # 3rd draw
        loc = kMC.next_position_draw(lat, reac, ersl, occ, r)
        
        lat, occ, ersl = kMC.update_lattice_configuration(lat, occ, ersl, reac, loc, DIRECTIONALITY)
        tick += 1
        
        """ Monitor the simulation progression every 500 loops
        """
        if tick%500 == 0: 
            compu_t = time.time()-t_init
            out = open('lattice_'+'%.f'%BOX_LENGTH+'_set_'+'%.f'%ite+'_output.txt','w')
            out.write('loop='+'%.f'%tick+', ite='+'%.2f'%((i/N_stat)*100)+', t_simu ='+'%.2f'%(compu_t/60)+' mins \n')
            out.close()
        
    t_comput = time.time() - t_init # time to compute the script
    data[i, 0] = t_comput/tick # velocity loop 
    data[i, 1] = t/tick # mean time step 

np.savetxt('lattice_'+'%.f'%BOX_LENGTH+'_set_'+'%.f'%ite+"_efficiency_measurements.dat", data)
