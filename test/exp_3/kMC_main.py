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
rd.seed(108) # Seed initialization 
tick_max = 10**4 # Stop condiontion on loops 
time_step = np.zeros(tick_max)
""" Parameters for initial lattice configuration 
     and for the motor walk model 
"""
t = 0 # Initial time
n = 125*10 # Nbr of dimer 
m = 13 # Nbr of Protofilaments 

Nr = 3 # Nbr of rate constants
kw = 100 # Cte rate of walking process  
km = 1 # Cte rate of detachment process 
rho = 0.1
kp = tl.attachment_rate_cte(rho, km, 0) # Cte rate of attachment process 

""" Initialize the lattice + event lists 
    lat: lattice of size (n,m)
    k: rate constant vector
    occ: number of enabled site for each reaction type, it's the occurence vector
    ersl[j] = list of the enable site positions for reaction j 
    lsre[k] = inverse of ersl = list of the reaction enable for site with index k 
    
"""
lat = kMC.latice_init(n, m)
k = kMC.rate_constant_init(Nr, kp, kw, km)
occ, ersl, lsre = kMC.initialize_eventlist(lat, Nr)

""" Initialized observables 
"""
tick = 0 # number of loops 
tick_visu = 250 # number of loops between 2 frames
bar = progressbar.ProgressBar(max_value=tick_max, widgets=widgets).start()

""" Lattice visualisation over time
"""
# fig, ax = plt.subplots(figsize=(12,12)) # Visu
# title = ax.text(0.5,1.25, "", bbox={'facecolor':'w','edgecolor':'black', 'alpha':0.75, 'pad':5},
#                 transform=ax.transAxes, ha="center")
# title.set_text(r"# Loops = {},".format(tick)+" t = {}s".format(np.round(t,2)))
# im = ax.imshow(np.swapaxes(lat,1,0), interpolation='nearest',cmap='tab20c_r',animated=True,vmin=-6,vmax=6)
# ax.set_ylim(0,.5)
# ax.set_xlim(0,n)
# plt.axis('on')
# plt.draw()
# plt.tight_layout()

""" loop the simulation while the stop condition is not verified
"""
while tick < tick_max: 
    r = rd.random() # 1st draw 
    C, reac = kMC.next_reaction_draw(k, occ, r)
    
    r = rd.random() # 2nd draw
    t += -(np.log(r)/C) # Simulated time evolution
    time_step[tick] = -(np.log(r)/C)
    r = rd.randint(0,int(occ[reac])-1) # 3rd draw
    loc = ersl[reac][r]
    
    lat, occ, ersl, lsre = kMC.update_lattice_configuration(lat, occ, ersl, lsre, reac, loc)
    tick += 1
    
    """ lattice visualisation update
    """
    # if tick%tick_visu == 0:
    #     plt.cla()
    #     title = ax.text(0.5, 3, "", bbox={'facecolor':'w','edgecolor':'black', 'alpha':0.75, 'pad':5},
    #             transform=ax.transAxes, ha="center")
    #     title.set_text(r"# Loops= {},".format(tick)+" t = {} s".format(np.round(t,2)))
    #     ax.imshow(np.swapaxes(lat,1,0), interpolation='nearest',cmap='tab20c',animated=True,vmin=-6,vmax=6)
    #     plt.axis('on')
    #     plt.draw()
    #     plt.pause(.1)
    bar.update(tick)

t_comput = time.time() - t_init # time to compute the script
v_loop = t_comput/tick
print()
print("it takes", round(v_loop*10**3, 3), "ms per loop")
print(time_step.mean(), time_step.std())