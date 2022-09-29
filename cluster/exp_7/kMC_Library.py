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


def initial_event_list_reaction(i, j, lattice, occ, ersl, lsre, DIRECTIONALITY):
    """
    test which reactions is enabled at site at (i,j) 
    and update occ & ersl & lsre for the initial lattice configuration
    """
    L = len(lattice[0,:])
    l = len(lattice[:,0])
    x = int(j+i*L)
    m = DIRECTIONALITY
        
    if 3 < i < l-3: 
        if lattice[i,j] == 0:
            """ it's a defect """
            ngbr = tl.ngbr_s3(i, j, lattice) 
            a = ngbr[0] ; at = ngbr[1]
            b = ngbr[2] ; bt = ngbr[3]  
            N_ngbr = a + at + (b+bt)/2 
            if (0 < N_ngbr < 4 or bt > 0) and np.abs(lattice[i-m,j]) != 3 and np.abs(lattice[i-m,j]) != 4: 
                """ any defect with at least one non-vacant neighbor 
                    without a motor upstream can be incorporated a GTP-tubulin
                """
                occ[5] += 1
                ersl[5].append([i,j])
                lsre[x].append(5)
        
        if lattice[i,j] != 0 and lattice[i,j]%2 == 0: 
            """ any excited dimer can relax to its fondamental level
            """
            occ[6] += 1
            ersl[6].append([i,j])
            lsre[x].append(6)
                
        if lattice[i,j] < 0 and lattice[i+1,j] != 0: 
            """ any GTP-dimer incorporated into the lattice can be hydrolized 
                incorporated = presence of 1 neighbor in the (+) direction
            """
            occ[4] += 1
            ersl[4].append([i,j])
            lsre[x].append(4)

        if  0 < np.abs(lattice[i,j]) <= 2: 
            """ it's a dimer without a motor can detach from the lattice 
            """
            ngbr = tl.ngbr_s3(i, j, lattice) 
            a = ngbr[0] ; at = ngbr[1]
            b = ngbr[2] ; bt = ngbr[3]  
            ind = int(a*(5*5*3) + at*(5*5) + b*5 + bt) 
            if lattice[i,j]%2 != 0: 
                # non-excited state 
                occ[7+ind] += 1
                ersl[7+ind].append([i,j])
                lsre[x].append(7+ind)
            else: 
                # excited state 
                occ[7+225+ind] += 1
                ersl[7+225+ind].append([i,j])
                lsre[x].append(7+225+ind)
            if 0 < np.abs(lattice[i+m,j]) <= 2: 
                """ when the adjacent dimer is without motor, a new motor can attach
                """
                occ[0] += 1
                ersl[0].append([i,j])
                lsre[x].append(0)

        if 2 < np.abs(lattice[i,j]) <= 4: 
            """ it's a dimer with a rear-head above 
            """
            if 4 < abs(lattice[i+m,j]) <= 6: 
                """ when the 2-head are attached, the motor can detach 
                """
                occ[1] += 1
                ersl[1].append([i,j])
                lsre[x].append(1)
                if 0 <= np.abs(lattice[i+2*m,j]) <= 2: 
                    """ when the dimer next to the front head is vacant, 
                    the motor can walk 
                    """
                    occ[3] += 3
                    ersl[3].append([i,j])
                    lsre[x].append(3)
            elif lattice[i+m,j] == 0: 
                """ if only one head is bound to the lattice, then the 
                    enhanced detachment is enabled and the below dimer can 
                    detach according to its state
                """
                occ[2] += 2
                ersl[2].append([i,j])
                lsre[x].append(2)
                
                ngbr = tl.ngbr_s3(i, j, lattice) 
                a = ngbr[0] ; at = ngbr[1]
                b = ngbr[2] ; bt = ngbr[3]  
                ind = int(a*(5*5*3) + at*(5*5) + b*5 + bt) 
                
                if lattice[i,j]%2 != 0: 
                    # non-excited dimer
                    occ[7+ind] += 1
                    ersl[7+ind].append([i,j])
                    lsre[x].append(7+ind)
                else: 
                    # excited state 
                    occ[7+225+ind] += 1
                    ersl[7+225+ind].append([i,j])
                    lsre[x].append(7+225+ind)        
                
    return (occ, ersl, lsre)


def initialize_eventlist(lat, Nbr_reac, DIRECTIONALITY): 
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
            occ, ersl, lsre = initial_event_list_reaction(i, j, lat, occ, ersl, lsre, DIRECTIONALITY)
    return (occ, ersl, lsre)



def update_enabled_reaction(i, j, lattice, DIRECTIONALITY):
    """ test which reactions are enabled for the site (i, j)
    """
    lsre_x = []  
    L = len(lattice[0,:])
    l = len(lattice[:,0])
    m = DIRECTIONALITY
        
    if 3 < i < l-3: 
        if lattice[i,j] == 0:
            """ it's a defect """
            ngbr = tl.ngbr_s3(i, j, lattice) 
            a = ngbr[0] ; at = ngbr[1]
            b = ngbr[2] ; bt = ngbr[3]  
            N_ngbr = a + at + (b+bt)/2 
            
            if (0 < N_ngbr < 4 or bt > 0) and np.abs(lattice[i-m,j]) != 3 and np.abs(lattice[i-m,j]) != 4: 
                """ any defect with at least one non-vacant neighbor 
                    without a motor upstream can be incorporated a GTP-tubulin
                """
                lsre_x.append(5)
        
        if lattice[i,j] != 0 and lattice[i,j]%2 == 0: 
            """ any excited dimer can relax to its fondamental level
            """
            lsre_x.append(6)
                
        if lattice[i,j] < 0 and lattice[i+1,j] != 0: 
            """ any GTP-dimer incorporated into the lattice can be hydrolized 
                incorporated = presence of 1 neighbor in the (+) direction
            """
            lsre_x.append(4)

        if  0 < np.abs(lattice[i,j]) <= 2: 
            """ it's a dimer without a motor can detach from the lattice 
            """
            ngbr = tl.ngbr_s3(i, j, lattice) 
            a = ngbr[0] ; at = ngbr[1]
            b = ngbr[2] ; bt = ngbr[3]      
            ind = int(a*(5*5*3) + at*(5*5) + b*5 + bt) 
            if lattice[i,j]%2 != 0: 
                # non-excited state 
                lsre_x.append(7+ind)
            else: 
                # excited state 
                lsre_x.append(7+225+ind)
            if 0 < np.abs(lattice[i+m,j]) <= 2: 
                """ when the adjacent dimer is without motor, a new motor can attach
                """
                lsre_x.append(0)

        if 2 < np.abs(lattice[i,j]) <= 4: 
            """ it's a dimer with a rear-head above 
            """
            if 4 < abs(lattice[i+m,j]) <= 6: 
                """ when the 2-head are attached, the motor can detach 
                """
                lsre_x.append(1)
                if 0 <= np.abs(lattice[i+2*m,j]) <= 2: 
                    """ when the dimer next to the front head is vacant, 
                    the motor can walk 
                    """
                    lsre_x.append(3)
            elif lattice[i+m,j] == 0: 
                """ if only one head is bound to the lattice, then the 
                    enhanced detachment is enabled and the below dimer can 
                    detach according to its state
                """
                lsre_x.append(2)
                
                ngbr = tl.ngbr_s3(i, j, lattice) 
                a = ngbr[0] ; at = ngbr[1]
                b = ngbr[2] ; bt = ngbr[3]     
                ind = int(a*(5*5*3) + at*(5*5) + b*5 + bt) 
                
                if lattice[i,j]%2 != 0: 
                    # non-excited dimer
                    lsre_x.append(7+ind)
                else: 
                    # excited state 
                    lsre_x.append(7+225+ind)    

    return lsre_x
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def scan(lsre, ersl, occ, loc, lat, DIRECTIONALITY):
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
            lsre_x = update_enabled_reaction(I, J, lat, DIRECTIONALITY)
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


def update_lattice_configuration(lat, occ, ersl, lsre, reac, loc, DIRECTIONALITY):
    """ change the lattice configuration with new reaction 
        and create corresponding event lists (ersl + occ + lsre)
    """
    x = int(loc[0])
    y = int(loc[1])
    m = DIRECTIONALITY
    
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
    occ, ersl, lsre = scan(lsre, ersl, occ, loc, lat, DIRECTIONALITY)    
    return (lat, occ, ersl, lsre)


