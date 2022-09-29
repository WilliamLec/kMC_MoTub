"""
kMC functions for the tradtionnal VSS method 
for a simple model of 2-head motors walk 
"""
import numpy as np
import Tools as tl

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
    #lattice[DEFECT_X_POSITION, DEFECT_Y_POSITION] = 0
    
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

def reaction_test(i, j, lattice, occ, ersl, TABLE_NO_DEFECT, TABLE_DEFECT, D):
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
                ersl[x, 1] += 1
                
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
                ersl[x, 2+ind] += 1

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
                ersl[x, 2+60+ind] += 1
                
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
                ersl[x, 2+60+30+ind] += 1
        
        if lattice[i, j] < 0 and lattice[i+1, j] != 0:
            """ it's a GTP-tubulin incorporated within the lattice """
            occ[0] += 1
            ersl[x, 0] += 1

    return occ, ersl


def scan(lat, Nr, DIRECTIONALITY, TABLE_NO_DEFECT, TABLE_DEFECT): 
    """ scan the whole lattice to obtain the enable reaction site list (ersl) 
            + occurence vector (occ)
    """
    
    n = len(lat[:,0]) 
    m = len(lat[0,:])
    D = DIRECTIONALITY
    ersl = np.zeros((n*m, Nr)) 
    occ = np.zeros(Nr) 

    for i in range(n):
        for j in range(m):
            occ, ersl = reaction_test(i, j, lat, occ, ersl, TABLE_NO_DEFECT, TABLE_DEFECT, D)
    return occ, ersl


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
    return C, reac


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
    

def update_lattice_configuration(lat, occ, ersl, reac, loc, DIRECTIONALITY, TABLE_DEFECT, TABLE_NO_DEFECT):
    """ change the lattice configuration with new reaction 
        and create corresponding event lists (ersl + occ)
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
    Nr = len(occ)
    occ, ersl = scan(lat, Nr, DIRECTIONALITY, TABLE_NO_DEFECT, TABLE_DEFECT)     
    return lat, occ, ersl
