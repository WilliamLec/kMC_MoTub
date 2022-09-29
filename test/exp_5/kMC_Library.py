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


def rate_constant_initizalition(k_MOTOR_ATTACH, k_MOTOR_DETACH, 
                                ENHANCED_FACTOR, k_WALK, 
                                HYDROLYSIS_TIME, k_DIMER_ATTACH, 
                                RELAXATION_TIME, 
                                LONGITUDINAL_ENERGY_GDP_TUBULIN,
                                LATERAL_ENERGY_GDP_TUBULIN,
                                LONGITUDINAL_ENERGY_GTP_TUBULIN,
                                LATERAL_ENERGY_GTP_TUBULIN, 
                                MOTOR_DESTABILIZATION):
    """ Return the rate constant vector for MoTub model
        
        There is 458 different reaction rate values : 
        
        0. motor attachment 
        1. motor detachment
        2. enhanced motor detachment 
        3. motor walk 
        4. GTP-hydrolysis 
        5. GTP-tubulin attachment
        6. Dimer relaxation 
        7 up to 232. non-destabilized dimer detachment 
        233 up to 458. destabilized dimer detachment 
    """

    k = np.zeros(7)
    k[0] = k_MOTOR_ATTACH 
    k[1] = k_MOTOR_DETACH
    k[2] = k_MOTOR_DETACH*ENHANCED_FACTOR
    k[3] = k_WALK
    k[4] = HYDROLYSIS_TIME**-1
    k[5] = k_DIMER_ATTACH
    k[6] = RELAXATION_TIME**-1
    
    DGbt = LONGITUDINAL_ENERGY_GTP_TUBULIN + LATERAL_ENERGY_GTP_TUBULIN
    DGl_d = LONGITUDINAL_ENERGY_GDP_TUBULIN
    DGL_d = LATERAL_ENERGY_GDP_TUBULIN
    DGl_t = LONGITUDINAL_ENERGY_GTP_TUBULIN
    DGL_t = LATERAL_ENERGY_GTP_TUBULIN
    DGp = MOTOR_DESTABILIZATION
    
    k_DETACH = np.zeros(3*3*5*5)
    k_DETACH_DESTABILIZED = np.zeros(3*3*5*5)
    index = 0 
    for a in range(3):
        for at in range(3):
            for b in range(5):
                for bt in range(5):   
                    DGb = a*DGl_d + at*DGl_t + (b/2)*DGL_d + (bt/2)*DGL_t
                    k_DETACH[index] = np.exp(DGbt - DGb)
                    k_DETACH_DESTABILIZED[index] = np.exp(DGbt - DGb + DGp)
                    index +=1 
    k = np.append(k, k_DETACH, axis=0) 
    k = np.append(k, k_DETACH_DESTABILIZED, axis=0)
    return k 

def reac(i, j, lattice, occ, ersl, m):
    """
    test which reactions is enabled at site at (i,j) 
    and return the updated occ & ersl
    """
    l = len(lattice[:,0])
    L = len(lattice[0,:])
    x = int(j+i*L)
    ##-------------------------------------------------------------------------
    if 3 < i < l-3: 
        if lattice[i,j] == 0:
            """ it's a defect """
            ngbr = tl.ngbr_s3(i, j, lattice) ; N_ngbr = sum(ngbr)
            a = ngbr[0] ; at = ngbr[1]
            b = ngbr[2] ; bt = ngbr[3]
            
            if (0 < N_ngbr < 4 or bt > 0) and np.abs(lattice[i-m,j]) != 3 and np.abs(lattice[i-m,j]) != 4: 
                """ any defect with at least one non-vacant neighbor 
                    without a motor upstream can be incorporated a GTP-tubulin
                """
                occ[5] += 1
                ersl[x, 5] += 1
                
        if lattice[i,j] != 0 and lattice[i,j]%2 == 0: 
            """ any excited dimer can relax to its fondamental level
            """
            occ[6] += 1
            ersl[x, 6] += 1
            
        if lattice[i,j] < 0 and lattice[i+1,j] != 0: 
            """ any GTP-dimer incorporated into the lattice can be hydrolized 
                incorporated = presence of 1 neighbor in the (+) direction
            """
            occ[4] += 1
            ersl[x, 4] += 1
                
        if  0 < np.abs(lattice[i,j]) <= 2: 
            """ it's a dimer without a motor can detach from the lattice 
            """
            ngbr = tl.ngbr_s3(i, j, lattice) ; N_ngbr = sum(ngbr)
            a = ngbr[0] ; at = ngbr[1]
            b = ngbr[2] ; bt = ngbr[3]      
            ind = int(a*(5*5*3) + at*(5*5) + b*5 + bt) 
            if lattice[i,j]%2 != 0: 
                # non-excited state 
                occ[7+ind] += 1
                ersl[x, 7+ind] += 1 
            else: 
                # excited state 
                occ[7+225+ind] += 1
                ersl[x, 7+225+ind] += 1 
                
            if 0 < np.abs(lattice[i+m,j]) <= 2: 
                """ when the adjacent dimer is without motor, a new motor can attach
                """
                occ[0] += 1
                ersl[x, 0] += 1
        
        if 2 < np.abs(lattice[i,j]) <= 4: 
            """ it's a dimer with a rear-head above 
            """
            if 4 < abs(lattice[i+m,j]) <= 6: 
                """ when the 2-head are attached, the motor can detach 
                """
                occ[1] += 1
                ersl[x, 1] += 1
                
                if 0 <= np.abs(lattice[i+2*m,j]) <= 2: 
                    """ when the dimer next to the front head is vacant, 
                    the motor can walk 
                    """
                    occ[3] += 1
                    ersl[x, 3] += 1
                    
            elif lattice[i+m,j] == 0: 
                """ if only one head is bound to the lattice, then the 
                    enhanced detachment is enabled and the below dimer can 
                    detach according to its state
                """
                occ[2] += 1
                ersl[x, 2] += 1
                
                ngbr = tl.ngbr_s3(i, j, lattice) ; N_ngbr = sum(ngbr)
                a = ngbr[0] ; at = ngbr[1]
                b = ngbr[2] ; bt = ngbr[3]     
                ind = int(a*(5*5*3) + at*(5*5) + b*5 + bt) 
                
                if lattice[i,j]%2 != 0: 
                    # non-excited dimer
                    occ[7+ind] += 1
                    ersl[x, 7+ind] += 1 
                else: 
                    # excited tubulin 
                    occ[7+225+ind] += 1
                    ersl[x, 7+225+ind] += 1 
    return occ, ersl


def scan(lat, Nr, DIRECTIONALITY): 
    """ scan the whole lattice to obtain the enable reaction site list (ersl) 
            + occurence vector (occ)
    """
    
    n = len(lat[:,0]) 
    m = len(lat[0,:])
    ersl = np.zeros((n*m, Nr)) 
    occ = np.zeros(Nr) 

    for i in range(n):
        for j in range(m):
            occ, ersl = reac(i, j, lat, occ, ersl, DIRECTIONALITY)
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
    

def update_lattice_configuration(lat, occ, ersl, reac, loc, m):
    """ change the lattice configuration with new reaction 
        and create corresponding event lists (ersl + occ)
    """
    x = int(loc[0])
    y = int(loc[1])
    
    """ check if there is some defect around
    """
    if lat[x,y] != 0: 
        s = abs(lat[x,y])/lat[x,y]
    if lat[x+m,y] != 0: # not a defect
        s1 = abs(lat[x+m,y])/lat[x+m,y]
    if lat[x+2*m,y] != 0: # not a defect
        s2 = abs(lat[x+2*m,y])/lat[x+2*m,y]
        
        
    if reac == 0:  
        """ a motor attaches the lattice 
        """
        lat[x,y] += 2*s 
        lat[x+m,y] += 4*s1 
        
    elif reac == 1: 
        """ a motor with 2-heads bound detaches the lattice
        """
        lat[x,y] -= 2*s 
        lat[x+m,y] -= 4*s1 
        
    elif reac == 2: 
        """ a motor with 1-head bound detaches the lattice
        """
        lat[x,y] -= 2*s 
    
    elif reac == 3: 
        """ a motor take a step 
        """
        if lat[x+2*m,y] != 0: 
            """ there is no defect ahead
            """
            lat[x,y] -= 2*s 
            lat[x+m,y] -= 2*s1 
            if lat[x+m,y]%2 == 1: 
                """ the dimer below the front head get destabilized if not already
                """
                lat[x+m,y] += s1
            lat[x+2*m,y] += 4*s2 
        else: 
            """ the motor is next to the MT edge or a defect, so only 1 head 
                is let bound to the lattice
            """
            lat[x,y] -= 2*s 
            lat[x+m,y] -= 2*s1 
            if lat[x+m,y]%2 == 1: 
                """ the dimer below the front head get destabilized if not already
                """
                lat[x+m,y] += s1            
    elif reac == 4: 
        """ The GTP-dimer get hydrolized 
        """
        lat[x,y] = np.abs(lat[x,y])
    
    elif reac == 5: 
        """ a GTP-tubulin attaches to the lattice
        """
        lat[x,y] = -1
    
    elif reac == 6: 
        """ the excited dimer relaxes 
        """
        lat[x,y] -= s
    
    elif reac >= 7: 
        """ the dimer detaches from the lattice 
        """
        lat[x,y] = 0      

    """ create the new ersl + occ 
    """
    Nr = len(occ)
    occ, ersl = scan(lat, Nr, m)         
    return lat, occ, ersl
