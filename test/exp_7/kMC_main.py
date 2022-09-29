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
import matplotlib.pyplot as plt
import kMC_Library as kMC
import Tools as tl
import progressbar
  
widgets = [' [',
          progressbar.Timer(format= 'elapsed time: %(elapsed)s'),
          '] ',
            progressbar.Bar('*'),' (',
            progressbar.ETA(), ') ',
          ]

def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

t_init = time.time()
rd.seed(18) # Seed initialization 
tick_max = 1*10**4 # Stop condiontion on loops 

""" Parameters for initial lattice configuration 
     and for the MoTub model 
"""
t = 0 # Initial time
BOX_LENGTH = 200 
BOX_WIDTH = 13 
LATTICE_HEIGHT = 125# int(0.5*BOX_LENGTH)
SEED_HEIGHT = 3
CAP_HEIGHT = 3

DEFECT_X_POSITION = int(0.5*LATTICE_HEIGHT)
DEFECT_Y_POSITION = int(0.5*BOX_WIDTH)

Nr = 7 + (3*3*5*5)*2 # Nbr of rate constants
k_MOTOR_DETACH = 1 # in s-1 
ENHANCED_FACTOR = 100  
MOTOR_DENSITY = 0. # motor density
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

MOTOR_DESTABILIZATION = 2 # in kT
DIRECTIONALITY = 1 

""" Initialize the lattice + event lists 
    lat: lattice of size (n,m)
    k: rate constant vector
    occ: number of enabled site for each reaction type, it's the occurence vector
    ersl: enabled reaction site list, it's an one-hot encoding matrice     
"""
lat = kMC.latice_initialization(BOX_LENGTH, BOX_WIDTH, LATTICE_HEIGHT, SEED_HEIGHT, CAP_HEIGHT,
                                DEFECT_X_POSITION, DEFECT_Y_POSITION)
k = kMC.rate_constant_initizalition(k_MOTOR_ATTACH, k_MOTOR_DETACH, ENHANCED_FACTOR, k_WALK,
                                    HYDROLYSIS_TIME, k_DIMER_ATTACH, RELAXATION_TIME, LONGITUDINAL_ENERGY_GDP_TUBULIN,
                                    LATERAL_ENERGY_GDP_TUBULIN, LONGITUDINAL_ENERGY_GTP_TUBULIN,
                                    LATERAL_ENERGY_GTP_TUBULIN, MOTOR_DESTABILIZATION)
occ, ersl, lsre = kMC.initialize_eventlist(lat, Nr, DIRECTIONALITY)

""" Initialized observables 
"""
tick = 0 # number of loops 
tick_visu = 250 # number of loops between 2 frames
bar = progressbar.ProgressBar(max_value=tick_max, widgets=widgets).start()

""" Lattice visualisation over time
"""
fig, ax = plt.subplots(figsize=(12,12)) # Visu
title = ax.text(0.5,1.25, "", bbox={'facecolor':'w','edgecolor':'black', 'alpha':0.75, 'pad':5},
                transform=ax.transAxes, ha="center")
title.set_text(r"# Loops = {},".format(tick)+" t = {}s".format(np.round(t,2)))
im = ax.imshow(np.swapaxes(lat,1,0), interpolation='nearest',cmap='tab20c',animated=True,vmin=-1,vmax=1)
ax.set_ylim(0,.5)
ax.set_xlim(0,BOX_LENGTH)
plt.axis('on')
plt.draw()
plt.tight_layout()

""" loop the simulation while the stop condition is not verified
"""
while tick < tick_max: 
    r = rd.random() # 1st draw 
    C, reac = kMC.next_reaction_draw(k, occ, r)
    
    r = rd.random() # 2nd draw
    t += -(np.log(r)/C) # Simulated time evolution
    
    r = rd.randint(0,int(occ[reac])-1) # 3rd draw
    loc = ersl[reac][r]
    
    lat, occ, ersl, lsre = kMC.update_lattice_configuration(lat, occ, ersl, lsre, reac, loc, DIRECTIONALITY)
    tick += 1
    
    """ lattice visualisation update
    """
    if tick%tick_visu == 0:
        plt.cla()
        title = ax.text(0.5, 3, "", bbox={'facecolor':'w','edgecolor':'black', 'alpha':0.75, 'pad':5},
                transform=ax.transAxes, ha="center")
        title.set_text(r"# Loops= {},".format(tick)+" t = {} s".format(np.round(t,2)))
        ax.imshow(np.swapaxes(lat,1,0), interpolation='nearest',cmap='tab20c',animated=True,vmin=-6,vmax=6)
        plt.axis('on')
        plt.draw()
        plt.pause(.1)
    bar.update(tick)

t_comput = time.time() - t_init # time to compute the script
v_loop = t_comput/tick
print()
print("it takes", round(v_loop*10**3, 3), "ms per loop")
