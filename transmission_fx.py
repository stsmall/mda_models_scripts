# -*- coding: utf-8 -*-
import numpy as np

def initialize_dfHost(infhost, muTrans, sigma, sizeTrans):
    '''Creates a transmission matrix for locations of infected hosts

    Parameters
    ----------
    infhost: 
    muTrans:(float) 
        parameter of neg binomial
    sizeTrans:(float) 
        parameter of neg binomial
    sigma:(float) 
        parameter for calculating dispersal distance

    Returns
    -------
    transmission_matrix: dataframe
        polar coordinate positions for pops
    dispersal: float
    '''
    print("starting trans_init")
    #set dispersal parameters
    dispersal = 4*math.pi*sigma**2
    pos = np.random.negative_binomial(sizeTrans,
            sizeTrans/float((sizeTrans+muTrans), (infhost, 2)))
    # :TODO fast polar conversion
    pospolar = to_polar(pos)

    transmission_matrix = pd.DataFrame({
                                        'hostpop' : range(infhost),
                                        'angle': pospolar[:, 0], 
                                        'radius': pospolar[:, 1],},
                                       )
    #{'meta_1': {'pop_1': ['x,y',[0,1,1,1,1,1]],'pop_2': ['x,y',[0,1,1,1,0,0]]}}
    return transmission_matrix, dispersal
    
def transmission_fx(mpop, transmission_mat, meta_popdict, dispersal, L3trans, hostpopsize):
    '''transmission function for worms between hosts
    L3trans:(int) number of successful transmission events
    mpop:(int) current metapopulation
    '''
    trans = 0
    mpop = "meta_" + str(mpop+1)
    while trans < L3trans:
        dpop = random.choice(meta_popdict[mpop].keys()) #random host to donate MF
        while [len(meta_popdict[mpop][dpop][i]) is 0 for i in meta_popdict[mpop][dpop] if "MF" in i].count(True) is 12: # make ceratin the donating pop actually has MF to donate
            dpop = random.choice(meta_popdict[mpop].keys())
        if len(transmission_mat[mpop][dpop][1]) >= hostpopsize:
            prob_newinf = 0 #all hosts are infected so no new infections
        else:
            prob_newinf = (float(1)/len(transmission_mat[mpop][dpop][1])) #prob of new infection is determined by number of uninfected  hosts
        #new infection
        if random.uniform(0, 1) < prob_newinf:
            newpop = "pop_{}".format(len(meta_popdict[mpop].keys()) + 1) #add new infection as population + 1
            meta_popdict[mpop][newpop]["J_1"] = []  #initialize new entry
            rmf = random.randint(1, 12) #which age class of MF is the L3 from
            while len(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) is 0:
                rmf = random.randint(1, 12) # keep choosing until MF age pop is not empty
            dmf = random.choice(meta_popdict[mpop][dpop]["MF_{}".format(rmf)])
            meta_popdict[mpop][newpop]["J_1"].append(dmf)
            transmission_mat = new_infection(transmission_mat, mpop, dpop, newpop, dispersal) #find new transmission position
        #reinfection
        else:
            rpop = random.choice([i for i, x in enumerate(transmission_mat[mpop][dpop][1]) if x == 1]) + 1 #choose a random value that is 1 as the receiving pop
            if meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"] is not list:
                meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"] = []  #in the case that J_1 is empty
            rmf = random.randint(1, 12)
            while len(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) is 0: # in the case that MF is empty in the donating
                rmf = random.randint(1, 12)
            dmf = random.choice(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) #choose random mf
            meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"].append(dmf) #add dontated mf
        trans += 1
    return meta_popdict, transmission_mat

def new_infection(transmission_mat, village, dpop, newpop, dispersal):
    '''
    Parameters
    ----------
    transmission_mat: np.array
        pairwise eucileadean distance between people 
    village: int
        metapopulation index
    dpop: int
        donation population, where genotypes are coming
    newpop: strings
        receiving population
    this is run every time a new individual is infected to rebuild the transmission prob matrix, updates'''
    x1n = np.random.uniform(-1*math.sqrt(dispersal), math.sqrt(dispersal))
    x2n = np.random.uniform(-1*math.sqrt(dispersal), math.sqrt(dispersal))
    # :TODO transmission matrix to dataframe

    transmission_matrix = pd.DataFrame({'village' : [],
                                        'hostpop' : [],
                                        'angle': [], 
                                        'radius': [],},
                                        index = ['1_p1'])

    distances = transmission_matrix.ix[x, y].distance()




    transmission_mat[village][newpop] =\
    ["{},{}".format(int(transmission_mat[village][dpop][0].split(",")[0]) + x1n, int(transmission_mat[mpop][dpop][0].split(",")[1]) + x2n), []]
    for tpop in range(1, len(transmission_mat[mpop].keys())):
        dist = math.sqrt((float(transmission_mat[mpop]["pop_{}".format(tpop)][0].split(",")[0]) - x1n)**2 + (float(transmission_mat[mpop]["pop_{}".format(tpop)][0].split(",")[1]) - x2n)**2)
        if dist < dispersal:
            dist = 1
        else:
            dist = 0
        transmission_mat[mpop][newpop][1].append(dist)
        transmission_mat[mpop]["pop_{}".format(tpop)][1].append(dist)
    transmission_mat[mpop][newpop][1].append(1)
    return transmission_mat
    