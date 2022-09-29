"""
kMC functions for the tradtionnal VSS method 
for a simple model of 2-head motors walk 
"""
import numpy as np

def latice_init(n, m):
    """ initialize the first configuration of the 2D lattice of size (n,m)
        initial configuration: empty lattice 
    """
    lat = np.zeros((n, m)) 
    return lat


def rate_constant_init(Nr, kp, kw, km):
    """ initialize the rate constant vector of Nr reactions
        0. motor attachment 
        1. motor walking 
        2. motor detachment
    """
    k = np.zeros(Nr)
    k[0] = kp 
    k[1] = kw 
    k[2] = km 
    return k 

def reac(i, j, lat, occ, ersl):
    """
    test which reactions is enabled at site at (i,j) 
    and update occ & ersl
    """
    L = len(lat[0,:])
    l = len(lat[:,0])
    x = int(j+i*L) 
    if lat[i,j] == 0 and lat[(i+1)%l,j] == 0:
        """ 2 vacant sites that are adjacent -> motor can attach 
        """
        occ[0] += 1 
        ersl[x,0] = 1
    elif lat[i,j] == 1: 
        """ presence of rear head -> motor can detach 
        """
        occ[2] += 1 
        ersl[x,2] = 1 
        if lat[(i+2)%l,j] == 0:
            """ if site ahead is vacant -> motor can walk
            """
            occ[1] += 1 
            ersl[x, 1] = 1 
    return (occ,ersl)


def scan(lat, Nr): 
    """ scan the whole lattice to obtain the enable reaction site list (ersl) 
            + occurence vector (occ)
    """
    
    n = len(lat[:,0]) 
    m = len(lat[0,:])
    ersl = np.zeros((n*m, Nr)) 
    occ = np.zeros(Nr) 

    for i in range(n):
        for j in range(m):
            occ, ersl = reac(i, j, lat, occ, ersl)
    return (occ, ersl)


def next_reaction_draw(k, occ, r):
    """ compute a probability for each reaction, and select 1 at random
        return the normalization factor & next reaction type
    """
    Nr = len(k)      
    kn = np.multiply(k, occ)  
    P = np.zeros(Nr)
    for i in range(Nr):
        """ non-normalized cummulative probability 
        to define intervals for each reaction
        """
        P[i] = np.sum(kn[0:(i+1)%(Nr+1)]) 
    C = np.sum(kn) # Normalizing factor
    Pn = np.zeros(Nr+1) 
    for i in range(Nr):
        """ normalized cummulative probability interval
        """
        Pn[i+1]= P[i]/C
    for i in range(0,Nr):
        if Pn[i] <= r < Pn[i+1]: 
            """ next reaction will be reaction n°i
            """
            reac = i  
    return (C, reac)


def next_position_draw(lat, reac, ersl, occ, r):
    """ construct the list of all the site that enabled the next reaction
        by scanning the corresponding column of ersl and select 1 at random
        return the (x,y) position of next reaction
    """
    L = len(lat[0,:])  
    liste = np.zeros((int(occ[int(reac)]), 2))
    idx = 0
    for k in range(len(ersl[:,0])):
       if ersl[k,int(reac)] == 1:  
           liste[idx, 1] = k%(L)
           liste[idx, 0] = (k - k%L)/L  
           idx += 1 
    x = liste[r, 0] ; y = liste[r, 1]
    return (x,y)
    

def update_lattice_configuration(lat, occ, ersl, reac, loc):
    """ change the lattice configuration with new reaction 
        and create corresponding event lists (ersl + occ)
    """
    x = int(loc[0]) ; y = int(loc[1])
    l = len(lat[:, 0]) 
    
    if reac == 0:  
        """ a motor attaches the lattice 
        """
        lat[x, y] = 1  # rear head on the lattice 
        lat[(x+1)%l, y] = 2 # Leading head on the lattice
    elif reac == 1:
        """ a motor takes a step on the lattice
        toward the increasing x-axis 
        """
        lat[x, y] = 0
        lat[(x+1)%l, y] = 1  
        lat[(x+2)%l, y] = 2 
    elif reac == 2: 
        """ a motor detaches from the lattice
        """
        lat[x, y] = 0    
        lat[(x+1)%l, y] = 0        

    """ create the new ersl + occ 
    """
    Nr = len(occ)
    occ, ersl = scan(lat, Nr)         
    return (lat,occ,ersl)
