"""
kMC functions for VSS method + local scan + use of list of list
for a simple model (flux of processive 2-head motors)
"""
import numpy as np
import Tools as tl

""" magical function to create list of list
"""
def mklist(n):
    for _ in range(n):
        yield []


def latice_initialization(BOX_LENGTH, BOX_WIDTH, 
                LATTICE_HEIGHT, SEED_HEIGHT, CAP_HEIGHT,
                DEFECT_X_POSITION, DEFECT_Y_POSITION):
    """ Return the first configuration of a 2D lattice 
        
        it can either be a : 
            1. MT with an initial GMPCPP seed 
            2. MT stabilized at its end 
            3. MT with a certain length but not attached at the lattice base
        
        A lattice can start with an initial defect
        
        - State = 0 <=> empty site
        - State = 1 <=> tubulin-GDP
        - State = -1 <=> tubulin-GTP
    """
    
    # Case 1. 
    lattice = np.zeros((BOX_LENGTH, BOX_WIDTH))
    lattice[:LATTICE_HEIGHT, :] = 1 
    lattice[:SEED_HEIGHT, :] = -1 
    lattice[LATTICE_HEIGHT-CAP_HEIGHT:LATTICE_HEIGHT, :] = -1 

    # Case 2. 
    # lattice = np.ones((BOX_LENGTH, BOX_WIDTH))

    # Case 3. 
    # MINUS_TIP_POSITION = int(BOX_LENGTH*0.25) 
    # lattice = np.zeros((BOX_LENGTH, BOX_WIDTH))
    # lattice[MINUS_TIP_POSITION:LATTICE_HEIGHT, :] = 1 
    # lattice[MINUS_TIP_POSITION:MINUS_TIP_POSITION+SEED_HEIGHT, :] = -1 
    #lattice[LATTICE_HEIGHT-CAP_HEIGHT:LATTICE_HEIGHT, :] = -1 
    
    # There is a defect along the lattice
    # lattice[DEFECT_X_POSITION, DEFECT_Y_POSITION] = 0
    
    return lattice

def rate_constant_initizalition(HYDROLYSIS_TIME, k_DIMER_ATTACH, 
                                RELAXATION_TIME, 
                                LONGITUDINAL_ENERGY_GDP_TUBULIN,
                                LATERAL_ENERGY_GDP_TUBULIN,
                                LONGITUDINAL_ENERGY_GTP_TUBULIN,
                                LATERAL_ENERGY_GTP_TUBULIN, 
                                LAMBDA_NO_DEFECT, LAMBDA_DEFECT):
    """ Return the rate constant vector for MoTub model
        
        There is 124 different reaction rate values : 
        
        0. GTP-hydrolysis 
        1. GTP-tubulin attachment
        2. Dimer relaxation 
        3 up to 63. effective koff for dimer w/out defect on same pf
        63 up to 93 effective koff for dimer w/ defect upstream
        93 up to 123 effective koff for dimer w/ defect downstream
    """

    k = np.zeros(2)
    k[0] = HYDROLYSIS_TIME**-1
    k[1] = k_DIMER_ATTACH
    
    DGbt = LONGITUDINAL_ENERGY_GTP_TUBULIN + LATERAL_ENERGY_GTP_TUBULIN
    DGl_d = LONGITUDINAL_ENERGY_GDP_TUBULIN
    DGL_d = LATERAL_ENERGY_GDP_TUBULIN
    DGl_t = LONGITUDINAL_ENERGY_GTP_TUBULIN
    DGL_t = LATERAL_ENERGY_GTP_TUBULIN
    
    k_DETACH = np.zeros(60)
    index = 0 
    for ai,a in enumerate([0,0,1,2]):
        for b in range(0,5):
            for bt in range(0,4-b+1):
                if ai != 0:
                    DGb = a*DGl_d + (2-a)*DGl_t + (b/2)*DGL_d + (bt/2)*DGL_t
                else: 
                    DGb = (b/2)*DGL_d + (bt/2)*DGL_t
                k_DETACH[index] = np.exp(DGbt-DGb)*LAMBDA_NO_DEFECT
                index +=1 
    k = np.append(k,k_DETACH,axis=0)

    koff_up = np.zeros(30)
    index = 0       
    for a in reversed(range(0,2)):
        for b in range(0,5):
            for bt in range(0,4-b+1): 
                DGb = a*DGl_d + (1-a)*DGl_t + (b/2)*DGL_d + (bt/2)*DGL_t
                koff_up[index] = np.exp(DGbt-DGb)*LAMBDA_DEFECT
                index += 1
    k = np.append(k,koff_up,axis=0)
    
    koff_down = np.zeros(30)
    index = 0       
    for a in reversed(range(0,2)):
        for b in range(0,5):
            for bt in range(0,4-b+1): 
                DGb = a*DGl_d + (1-a)*DGl_t + (b/2)*DGL_d + (bt/2)*DGL_t
                koff_down[index] = np.exp(DGbt-DGb)
                index += 1
    k = np.append(k,koff_down,axis=0)
    
    return k 
                    
def initial_event_list_reaction(i, j, lattice, occ, ersl, lsre, TABLE_NO_DEFECT, TABLE_DEFECT, D):
    """
    test which reactions is enabled at site at (i,j) 
    and return the updated occ & ersl
    """
    l = len(lattice[:,0])
    L = len(lattice[0,:])
    x = int(j+i*L)
    ##-------------------------------------------------------------------------
    if 3 < i < l-3: 
        ngbr = tl.ngbr_s3(i, j, lattice) ; N_ngbr = sum(ngbr)
        a = ngbr[0] ; at = ngbr[1]
        b = ngbr[2] ; bt = ngbr[3]
        if lattice[i,j] == 0:
            """ it's a defect """
            if (0 < N_ngbr < 4 or bt > 0): 
                """ any defect with at least one non-vacant neighbor 
                    without a motor upstream can be incorporated a GTP-tubulin
                """
                occ[1] += 1
                ersl[1].append([i,j])
                lsre[x].append(1)
                
        if lattice[i,j] != 0: 
            """ it's a tubulin dimer
            """
            if a+at != 1: 
                """ it's a dimer far from defect or isolated dimer """
                found = False
                ind = len(TABLE_NO_DEFECT[:,0])-1
                while not(found): 
                    """ searching for the corresponding koff index"""
                    if np.all([a,at,b,bt] == TABLE_NO_DEFECT[ind,:]):
                        found = True
                    else: 
                        ind -= 1 
                occ[2+ind] += 1
                ersl[2+ind].append([i, j])
                lsre[x].append(2+ind)

            elif lattice[i+D, j] == 0: 
                """ there's a defect upstream """  
                found = False
                ind = len(TABLE_DEFECT[:,0])-1
                while not(found): 
                    """ searching for the corresponding koff index"""
                    if np.all([a,at,b,bt] == TABLE_DEFECT[ind,:]):
                        found = True
                    else: 
                        ind -= 1 
                occ[2+60+ind] += 1
                ersl[2+60+ind].append([i, j])
                lsre[x].append(2+60+ind)

                
            elif lattice[i-D, j] == 0: 
                """ there's a defect downstream """  
                found = False
                ind = len(TABLE_DEFECT[:,0])-1
                while not(found): 
                    """ searching for the corresponding koff index"""
                    if np.all([a,at,b,bt] == TABLE_DEFECT[ind,:]):
                        found = True
                    else: 
                        ind -= 1 
                occ[2+60+30+ind] += 1
                ersl[2+60+30+ind].append([i, j])
                lsre[x].append(2+60+30+ind)

        
        if lattice[i, j] < 0 and lattice[i+1, j] != 0:
            """ it's a GTP-tubulin incorporated within the lattice """
            occ[0] += 1
            ersl[0].append([i,j])
            lsre[x].append(0)
                       
    return (occ, ersl, lsre)


def initialize_eventlist(lat, Nbr_reac, TABLE_NO_DEFECT, TABLE_DEFECT, DIRECTIONALITY): 
    """ initialize the occurence array and the ersl list and lsre list
        occ[i] = # of site that enable the reaction i 
        ersl[j] = list of the enable site positions for reaction j 
        lsre[k] = inverse of ersl = list of the reaction enable for site with index k 
    """
    D = DIRECTIONALITY
    occ = np.zeros(Nbr_reac) 
    ersl = list(mklist(Nbr_reac))                                                     
    lsre = list(mklist(len(lat[:,0])*len(lat[0,:])))
    
    L = len(lat[0,:])
    l = len(lat[:,0])

    for i in range(l):
        for j in range(L):
            occ, ersl, lsre = initial_event_list_reaction(i, j, lat, occ, ersl, lsre, TABLE_NO_DEFECT, TABLE_DEFECT, D)
    return (occ, ersl, lsre)



def update_enabled_reaction(i, j, lattice, TABLE_NO_DEFECT, TABLE_DEFECT, DIRECTIONALITY):
    """ test which reactions are enabled for the site (i, j)
    """
    lsre_x = []  
    l = len(lattice[:,0])
    D = DIRECTIONALITY
     
    if 3 < i < l-3: 
        ngbr = tl.ngbr_s3(i, j, lattice) ; N_ngbr = sum(ngbr)
        a = ngbr[0] ; at = ngbr[1]
        b = ngbr[2] ; bt = ngbr[3]
        if lattice[i,j] == 0:
            """ it's a defect """
            if (0 < N_ngbr < 4 or bt > 0): 
                """ any defect with at least one non-vacant neighbor 
                    without a motor upstream can be incorporated a GTP-tubulin
                """
                lsre_x.append(1)
                
        if lattice[i,j] != 0: 
            """ it's a tubulin dimer
            """
            if a+at != 1: 
                """ it's a dimer far from defect or isolated dimer """
                found = False
                ind = len(TABLE_NO_DEFECT[:,0])-1
                while not(found): 
                    """ searching for the corresponding koff index"""
                    if np.all([a,at,b,bt] == TABLE_NO_DEFECT[ind,:]):
                        found = True
                    else: 
                        ind -= 1 
                lsre_x.append(2+ind)

            elif lattice[i+D, j] == 0: 
                """ there's a defect upstream """  
                found = False
                ind = len(TABLE_DEFECT[:,0])-1
                while not(found): 
                    """ searching for the corresponding koff index"""
                    if np.all([a,at,b,bt] == TABLE_DEFECT[ind,:]):
                        found = True
                    else: 
                        ind -= 1 
                lsre_x.append(2+60+ind)

                
            elif lattice[i-D, j] == 0: 
                """ there's a defect downstream """  
                found = False
                ind = len(TABLE_DEFECT[:,0])-1
                while not(found): 
                    """ searching for the corresponding koff index"""
                    if np.all([a,at,b,bt] == TABLE_DEFECT[ind,:]):
                        found = True
                    else: 
                        ind -= 1 
                lsre_x.append(2+60+30+ind)

        if lattice[i, j] < 0 and lattice[i+1, j] != 0:
            """ it's a GTP-tubulin incorporated within the lattice """
            lsre_x.append(0)   

    return lsre_x
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def scan(lsre, ersl, occ, loc, lat, TABLE_NO_DEFECT, TABLE_DEFECT, DIRECTIONALITY):
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
    half_width = 1
    
    for I in range(i-half_length,i+half_length+1):
        for J in range(j-half_width,j+half_width+1):
            I = I%l ; J = J%L
            lsre_x = update_enabled_reaction(I, J, lat, TABLE_NO_DEFECT, TABLE_DEFECT, DIRECTIONALITY)
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


def update_lattice_configuration(lat, occ, ersl, lsre, reac, loc, TABLE_NO_DEFECT, TABLE_DEFECT, DIRECTIONALITY):
    """ change the lattice configuration with new reaction 
        and create corresponding event lists (ersl + occ + lsre)
    """
    x = int(loc[0])
    y = int(loc[1])
             
    if reac == 0: 
        """ The GTP-dimer get hydrolized 
        """
        lat[x,y] = 1
    
    elif reac == 1: 
        """ a GTP-tubulin attaches to the lattice
        """
        lat[x,y] = -1
    
    elif reac > 1: 
        """ the dimer detaches from the lattice 
        """
        lat[x,y] = 0      

    """ create the new ersl + occ 
    """
    occ, ersl, lsre = scan(lsre, ersl, occ, loc, lat, TABLE_NO_DEFECT, TABLE_DEFECT, DIRECTIONALITY)
    return (lat, occ, ersl, lsre)


