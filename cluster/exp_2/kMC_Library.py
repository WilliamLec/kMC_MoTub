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

def reac(i, j, lat, occ, ersl, is_scan, idx):
    """
    test which reactions is enabled at site at (i,j) 
    and update occ & ersl
    """
    L = len(lat[0,:])
    l = len(lat[:,0])
    
    """ determine if it is used by the scan 
    """
    if not(is_scan):
        x = int(j+i*L)
    else:
        x = int(idx)
        
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

def initialize_event_lists(lat, Nr): 
    """ initialize the occurence array and the ersl matrice
    """
    
    n = len(lat[:,0]) 
    m = len(lat[0,:])
    ersl = np.zeros((n*m, Nr)) 
    occ = np.zeros(Nr) 

    for i in range(n):
        for j in range(m):
            occ, ersl = reac(i, j, lat, occ, ersl, False, 0)
    return (occ, ersl)


def scan(lat, loc, ersl, occ): 
    """ scan the neighborhood around
    the location of the new reaction site 
    to update the enable reaction site list (ersl) 
            + occurence vector (occ)
    the size of the scan region is determined by its half_width and half_length
    as ( 2 half_width + 1 x 2 half_lenght + 1)
    """
    
    x = int(loc[0])     # Reac x-position 
    y = int(loc[1])     # Reac y-position 
    l = len(lat[:,0]) # Nbr of dimer per pf 
    L = len(lat[0,:]) # Nbr of pf per microtubule 
    half_width = 1
    half_length = 2 
    if L == 1: 
        nbr_neighbors = (2*half_length + 1)
        half_width = 0
    elif L == 2: 
        nbr_neighbors = 2*(2*half_length + 1)
        half_width = 1
    elif L > 2: 
        nbr_neighbors = (2*half_width + 1)*(2*half_length + 1)
        
    """ create the ersl of the intial ersl of the neighborhood before 
    the reaction
    """
    ersl_neighborhood = np.zeros((nbr_neighbors,len(occ)))
    conv = np.zeros((nbr_neighbors,3))
    idx = 0 
    for a in range(-half_length,half_length+1):
        if L == 2:
            for b in range(2):
                conv[idx,0] = (x-a)%l
                conv[idx,1] = (y-b)%L
                conv[idx,2] = ((x-a)%l)*L+(y-b)%L
                ersl_neighborhood[idx,:] = ersl[int(conv[int(idx),2]),:]
                idx += 1 
        elif L != 2:
            for b in range(-half_width, half_width+1):
                conv[idx,0] = (x-a)%l
                conv[idx,1] = (y-b)%L
                conv[idx,2] = ((x-a)%l)*L+(y-b)%L
                ersl_neighborhood[idx,:] = ersl[int(conv[int(idx),2]),:]
                idx += 1 
    """substract from the occurence vector the occurence of the neighborhood
    """
    occ_neighborhood = np.zeros(len(occ))
    for idx1 in range(0,len(occ)):
        occ_neighborhood[idx1] = np.sum(ersl_neighborhood[:,idx1])
    occ -= occ_neighborhood 
    
    """ scan the neighborhood of the updated lattice 
    to create the occurence vector and ersl of the neighborhood 
    """
    occ_neighborhood = np.zeros(len(occ))
    ersl_neighborhood = np.zeros((15,len(occ)))
    for index in range(nbr_neighbors): 
        i = int(conv[index,0])
        j = int(conv[index,1])
        occ_neighborhood, ersl_neighborhood = reac(i, j, lat, occ_neighborhood, ersl_neighborhood, True, index)
    occ += occ_neighborhood
    for a in  range(nbr_neighbors):
        ersl[int(conv[a,2]),:] = ersl_neighborhood[a,:]     
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
    occ, ersl = scan(lat, loc, ersl, occ)     
    return (lat,occ,ersl)
