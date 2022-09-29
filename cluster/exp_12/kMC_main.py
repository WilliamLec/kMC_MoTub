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
t_init0 = time.time()

rd.seed(args.seed) # Seed initialization 
tick_max = 10**4 # Stop condiontion on loops 
N_stat = 10 
data = np.zeros((N_stat, 2))



""" Parameters for initial lattice configuration 
     and for the MoTub model 
"""
t = 0 # Initial time
BOX_LENGTH = args.L
BOX_WIDTH = 10
LATTICE_HEIGHT = int(0.5*BOX_LENGTH)
SEED_HEIGHT = 3
CAP_HEIGHT = 3

DEFECT_X_POSITION = int(0.5*LATTICE_HEIGHT)
DEFECT_Y_POSITION = int(0.5*BOX_WIDTH)

Nr = 2 + 60 + 30*2 # Nbr of rate constants
k_MOTOR_DETACH = 1 # in s-1 
ENHANCED_FACTOR = 100  
MOTOR_DENSITY = 0 # motor density
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
STABILIZATION_GTP_TUBULIN = 0.5*(CONSTRAST-1)*TOTAL_BINDING_ENERGY_GDP_TUBULIN
LONGITUDINAL_ENERGY_GTP_TUBULIN = LONGITUDINAL_ENERGY_GDP_TUBULIN + STABILIZATION_GTP_TUBULIN
LATERAL_ENERGY_GTP_TUBULIN = LATERAL_ENERGY_GDP_TUBULIN

""" Motor characteristics """ 
motor = "kinesin"
if  motor == "kinesin":
    DIRECTIONALITY = +1 # directionality of the motor, m = 1 (resp. -1) for + end (resp - end) 
elif motor == "dynein":
    DIRECTIONALITY = -1
rho = 0.0  # steady-density of motors
taur = .1 # relaxation time [s]
theta = 100 # omega'_d = theta x omega_d
DGp = 1 # motor penalty [kT]

""" Effective constants """ 
tabulation = np.loadtxt("./"+"%s"%motor+'_taur_'+'%.1f'%taur+'_theta_'+'%.f'%theta+'_tabulation_detachment_probabilities.dat')

if rho != 0: 
    cg_ND = tabulation[int(np.where(tabulation[:,0] == rho)[0]),1] 
    ce_ND = tabulation[int(np.where(tabulation[:,0] == rho)[0]),2]
    LAMBDA_NO_DEFECT = cg_ND + ce_ND*np.exp(DGp)  
    
    cg_D = tabulation[int(np.where(tabulation[:,0] == rho)[0]),3]
    ce_D = tabulation[int(np.where(tabulation[:,0] == rho)[0]),4]
    LAMBDA_DEFECT = cg_D + ce_D*np.exp(DGp)
else: 
    LAMBDA_DEFECT = 1
    LAMBDA_NO_DEFECT = 1

""" Table for koff index for tubulin-tubulin interaction configurations without defect """
TABLE_NO_DEFECT = np.zeros((60,4)) 
index = 0 
for ui,u in enumerate([0,0,1,2]):
    for v in range(0,5):
        for vt in range(0,4-v+1):  
            TABLE_NO_DEFECT[index,0] = u
            if ui != 0:
                TABLE_NO_DEFECT[index,1] = 2-u
            TABLE_NO_DEFECT[index,2] = v/2
            TABLE_NO_DEFECT[index,3] = vt/2
            index += 1
                
""" Table for koff index for tubulin-tubulin interaction configurations with longitudinal defect """
TABLE_DEFECT = np.zeros((30,4)) # all different configurations for tubulin-tubulin interaction with 1 missing longitudinal neighbor 
index = 0       
for u in reversed(range(0,2)):
    for v in range(0,5):
        for vt in range(0,4-v+1):  
            TABLE_DEFECT[index,0] = u
            TABLE_DEFECT[index,1] = 1-u
            TABLE_DEFECT[index,2] = v/2
            TABLE_DEFECT[index,3] = vt/2
            index += 1


""" loop over N_stat time for statistics purposes
"""
for i in range(N_stat):

    
    """ Initialize the lattice + event lists 
        lat: lattice of size (n,m)
        k: rate constant vector
        occ: number of enabled site for each reaction type, it's the occurence vector
        ersl: enabled reaction site list, it's an one-hot encoding matrice     
    """
    lat = kMC.latice_initialization(BOX_LENGTH, BOX_WIDTH, LATTICE_HEIGHT, SEED_HEIGHT, CAP_HEIGHT,
                                    DEFECT_X_POSITION, DEFECT_Y_POSITION)
    k = kMC.rate_constant_initizalition(HYDROLYSIS_TIME, k_DIMER_ATTACH, RELAXATION_TIME,
                                        LONGITUDINAL_ENERGY_GDP_TUBULIN, LATERAL_ENERGY_GDP_TUBULIN,
                                        LONGITUDINAL_ENERGY_GTP_TUBULIN, LATERAL_ENERGY_GTP_TUBULIN,
                                        LAMBDA_NO_DEFECT, LAMBDA_DEFECT)
    occ, ersl, lsre = kMC.initialize_eventlist(lat, Nr, TABLE_NO_DEFECT, TABLE_DEFECT, DIRECTIONALITY)
    

    """ Initialized observables 
    """
    t = 0 # Initial time
    t_init = time.time()
    tick = 0 # number of loops 
    
    
    """ loop the simulation while the stop condition is not verified
    """
    while tick < tick_max: 
        r = rd.random() # 1st draw 
        C, reac = kMC.next_reaction_draw(k, occ, r)
        
        r = rd.random() # 2nd draw
        t += -(np.log(r)/C) # Simulated time evolution
        
        r = rd.randint(0,int(occ[reac])-1) # 3rd draw
        loc = ersl[reac][r]
        
        lat, occ, ersl, lsre = kMC.update_lattice_configuration(lat, occ, ersl, lsre, reac, loc, TABLE_NO_DEFECT, TABLE_DEFECT, DIRECTIONALITY)
    
        """ Monitor the simulation progression every 500 loops
        """
        if tick%500 == 0: 
            compu_t = time.time()-t_init0
            out = open('lattice_'+'%.f'%BOX_LENGTH+'_set_'+'%.f'%ite+'_output.txt','w')
            out.write('loop='+'%.f'%tick+', ite='+'%.2f'%((i/N_stat)*100)+', t_simu ='+'%.2f'%(compu_t/60)+' mins \n')
            out.close()
        tick += 1
        
    t_comput = time.time() - t_init # time to compute the script
    data[i, 0] = t_comput/tick # velocity loop 
    data[i, 1] = t/tick # mean time step 

np.savetxt('lattice_'+'%.f'%BOX_LENGTH+'_set_'+'%.f'%ite+"_efficiency_measurements.dat", data)