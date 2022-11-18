"""
kMC functions for VSS method + local scan + use of list of list
for a simple model (flux of processive 2-head motors)
"""
import numpy as np

""" magical function to create list of list
"""
def mklist(n):
    for _ in range(n):
        yield []


def latice_init(n, m):
    """ initialize the first configuration of the 2D lattice of size (n,m)
        initial configuration: empty lattice 
    """
    
    lat = np.zeros((n, m)) 
    
    """ add a rectangle defect at mid length """
    
    ld = 5
    Ld = int(m/3)
    
    lat[int((n-ld)/2):int((n+ld)/2),int((m-Ld)/2):int((m+Ld)/2)] = -1 
    
    return lat

def rate_constant_init(Nr, kp, kw, km, theta):
    """ initialize the rate constant vector of Nr reactions
        0. motor attachment 
        1. motor walking 
        2. motor detachment
    """
    k = np.zeros(Nr)
    k[0] = kp 
    k[1] = kw 
    k[2] = km
    k[3] = km*theta
    return k 


def initial_event_list_reaction(i, j, lat, occ, ersl, lsre):
    """
    test which reactions is enabled at site at (i,j) 
    and update occ & ersl & lsre for the initial lattice configuration
    """
    L = len(lat[0,:])
    l = len(lat[:,0])
    x = int(j+i*L)

        
    if lat[i,j] == 0 and lat[(i+1)%l,j] == 0:
        """ 2 vacant sites that are adjacent -> motor can attach 
        """
        occ[0] += 1
        ersl[0].append([i,j])
        lsre[x].append(0)
        
    elif lat[i,j] == 1 and lat[(i+1)%l,j] == 2: 
        """ presence of both heads = "normal" detachment 
        """
        occ[2] += 1
        ersl[2].append([i,j])
        lsre[x].append(2)
        if lat[(i+2)%l,j] != 1:
            """ if site ahead is not a rear head motor -> motor can walk
                not rear head = either a vacant or a defect site
            """
            occ[1] += 1
            ersl[1].append([i,j])
            lsre[x].append(1)
        
    elif lat[i,j] == 1 and lat[(i+1)%l,j] == -1: 
        """ presence of only one head and defect ahead = "modified" detachment 
        """
        occ[3] += 1
        ersl[3].append([i,j])
        lsre[x].append(3)
        
    return (occ, ersl, lsre)


def initialize_eventlist(lat, Nbr_reac): 
    """ initialize the occurence array and the ersl list and lsre list
        occ[i] = # of site that enable the reaction i 
        ersl[j] = list of the enable site positions for reaction j 
        lsre[k] = inverse of ersl = list of the reaction enable for site with index k 
    """
    occ = np.zeros(Nbr_reac) 
    ersl = list(mklist(Nbr_reac))                                                     
    lsre = list(mklist(len(lat[:,0])*len(lat[0,:])))
    
    L = len(lat[0,:])
    l = len(lat[:,0])

    for i in range(l):
        for j in range(L):
            occ, ersl, lsre = initial_event_list_reaction(i, j, lat, occ, ersl, lsre)
    return (occ, ersl, lsre)



def update_enabled_reaction(i, j, lat):
    """ test which reactions are enabled for the site (i, j)
    """
    l = len(lat[:,0])
    lsre_x = []  
    
    if lat[i,j] == 0 and lat[(i+1)%l,j] == 0:
        """ 2 vacant sites that are adjacent -> motor can attach 
        """
        lsre_x.append(0)

    elif lat[i,j] == 1 and lat[(i+1)%l,j] == 2: 
        """ presence of both heads = "normal" detachment 
        """
        lsre_x.append(2)
        if lat[(i+2)%l,j] != 1:
            """ if site ahead is not a rear head motor -> motor can walk
                not rear head = either a vacant or a defect site
            """
            lsre_x.append(1)
        
    elif lat[i,j] == 1 and lat[(i+1)%l,j] == -1: 
        """ presence of only one head and defect ahead = "modified" detachment 
        """
        lsre_x.append(3)
    return lsre_x
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def scan(lsre, ersl, occ, loc, lat):
    """ scan the neighborhood around the location (i,j) of the new reaction site 
    to update the enable reaction site list (ersl) + occurence vector (occ)
    
    the size of the scan region is determined by its half_width and half_length
    as ( 2 half_width + 1 x 2 half_lenght + 1)
    """
    
    i = int(loc[0])    
    j = int(loc[1])    
    l = len(lat[:,0]) 
    L = len(lat[0,:]) 
    
    half_length = 3
    half_width = 0
    
    for I in range(i-half_length,i+half_length+1):
        for J in range(j-half_width,j+half_width+1):
            I = I%l ; J = J%L
            lsre_x = update_enabled_reaction(I, J, lat)
            lsre_x.sort()
            X = int(J+I*L)
            lsre[X].sort()
            
            """ check if the enabled reaction list at index x is different
            from the previous one
            
            if it is than update all the eventlist (occ, ersl, lsre) at this 
            index 
            """
            if lsre_x != lsre[X]: 
                """ add new enabled reaction     
                """
                for k in range(len(lsre_x)):
                    """ if reactions not already enabled then add them
                    """
                    if lsre[X].count(lsre_x[k]) == 0: # 
                        occ[lsre_x[k]] += 1
                        ersl[lsre_x[k]].append([I,J])
                        
                """ if reactions are not anymore enabled then revome them
                """
                for k in range(len(lsre[X])):
                    if lsre_x.count(lsre[X][k]) == 0: 
                        occ[lsre[X][k]] -= 1
                        ersl[lsre[X][k]].remove([I,J])
                        
                """ update the lsre list 
                """
                lsre[X] = lsre_x
    return (occ, ersl, lsre)


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


def update_lattice_configuration(lat, occ, ersl, lsre, reac, loc):
    """ change the lattice configuration with new reaction 
        and create corresponding event lists (lsre + ersl + occ)
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
        if lat[(x+2)%l, y] == -1: 
            """ presence of a defect ahead """
            lat[x, y] = 0
            lat[(x+1)%l, y] = 1  
        else : 
            lat[x, y] = 0
            lat[(x+1)%l, y] = 1  
            lat[(x+2)%l, y] = 2 
    
    elif reac == 2: 
        """ a motor detaches from the lattice
        """
        lat[x, y] = 0    
        lat[(x+1)%l, y] = 0   
        
    elif reac == 3: 
        """ one-headed motor detaches from the lattice
        """
        lat[x, y] = 0        

    """ update the event lists (lsre + ersl + occ)
    """
    occ, erls, lsre = scan(lsre, ersl, occ, loc, lat)
    return (lat, occ, ersl, lsre)

