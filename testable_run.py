# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import numpy as np
import pandas as pd
import math
from sklearn.metrics import pairwise_distances
import random
from scipy.stats import weibull_min
import subprocess

def agehost_fx(sex):
    '''calculate age of host and age of death
    Parameters
    ---------
    sex: str
         male or female

    Returns
    ------
    age: int
         age of host
    death: int
         age of death
    '''
    zero_10 = 0.24 #12%
    eleven_20 = 0.46 #11%
    twentyone_30 = 0.64 #9%
    thirtyone_40 = 0.78 #7%
    fourtyone_50 = 0.88 #5%
    fiftyone_60 = 0.94 #3%
    sixtyone_70 = 0.98 #2%                     
    #assign age
    agerand = np.random.random()
    if agerand <= zero_10:
         age = np.random.randint(1,10)
    elif agerand > zero_10 and agerand <= eleven_20:
         age = np.random.randint(11,20)
    elif agerand > eleven_20 and agerand <= twentyone_30:
         age = np.random.randint(21,30)
    elif agerand > twentyone_30 and agerand <= thirtyone_40:
         age = np.random.randint(31,40)
    elif agerand > thirtyone_40 and agerand <= fourtyone_50:
         age = np.random.randint(41,50)
    elif agerand > fourtyone_50 and agerand <= fiftyone_60:
         age = np.random.randint(51,60)
    elif agerand > fiftyone_60 and agerand <= sixtyone_70:
         age = np.random.randint(61,70)
    elif agerand > sixtyone_70:
         age = np.random.randint(71,80)
    
    #dictionary from actuarial tables
    
    deathdict = {}
    with open("act.tbl",'r') as tbl:
         for line in tbl:
              line = line.strip()
              deathdict["{}".format(line.split()[0])] = list(map(float,
                  line.split()[1:]))
    
    #when do they die
    death = deathdict[str(age)][int(sex)] + np.random.normal(0,6) + age
         
    return age, round(death)
    
def host_fx(villages, infhost, muTrans, sizeTrans):
    '''Creates a transmission matrix for locations of infected hosts
       default_test : host_fx(2, [100, 300], 100, 1)

    Parameters
    ----------
    infhost: list of pop. sizes
         number of infected hosts at intial time. prev*hostpopsize
    muTrans: float 
        parameter of neg binomial
    sizeTrans: float 
        parameter of neg binomial
    sigma: float
        dispersal distance in meters

    Returns
    -------
    dfHost: dataframe
    '''
    assert villages == len(infhost)
    coordinates = []
    host_idx = []
    for vill in range(villages):
         #list of host positions
         coordinates.extend(np.random.negative_binomial(sizeTrans, sizeTrans 
                              / float((sizeTrans+muTrans)), (infhost[vill],
                                  2)))
         for host in range(infhost[vill]):
             host_idx.append("v" + str(vill) + "h" + str(host + 1))
             
    sex = [random.choice("01") for i in range(sum(infhost))]
    age_death  = [agehost_fx(i) for i in sex]

    dfHost = pd.DataFrame({
                      'village' : np.repeat(range(villages), infhost),
                      'hostidx' : host_idx,
                      'sex' : sex,
                      'age' : [i[0] for i in age_death],
                      'agedeath' : [i[1] for i in age_death],
                      'coordinates' : coordinates,
                      'MDA' : np.zeros(sum(infhost)),
                      })
    dfHost = dfHost.loc[:, ['village', 'hostidx', 'sex',
            'age', 'agedeath', 'coordinates', 'MDA']]
              
    return dfHost
    
def coalsims_migmat_fx(villages, initial_migration, initial_distance_m, theta,
                     basepairs, mutation):
    '''Creates a string that represents a migration matrix between
    metapopulations. Migration here is a stepping stone not island model

    Format specified by ms (hudson 2000). Uses euclidian distances. The value
    of 2Nm is weighted by the distance from the next population as an
    exponential random variable.  The highest this can be is 2Nm/1

    Parameters
    ----------
    villages : int
        number of villages/metapopulations
    initial_migration : float
        migration for ms/scrm as the fraction of subpopulation i
        which is made up of migrants from subpopulation j each generation
    initial_distance_m : (list,float)
        distance between villages in meters
    theta : float, list
        theta values per village for locus
    basepairs : int
         length of locus
    mutation : float
         mutation rate of locus
         
    Returns
    -------
    Migration matrix: string
        migration matrix for scrm
    '''
    ne = theta / (4 * mutation * basepairs)

    if villages > 4:
        raise ValueError("only handles 4 villages ATM")
    elif villages < 4:
        if len(initial_distance_m) != ((villages) * (villages - 1) / 2):
            raise ValueError(("there are not adequate pairwise comparisons in"
                              "distance_m to match villages"))
        mig = []  # initiate blank migration list
        for meters in initial_distance_m:
            mig.append((initial_migration) / (np.random.exponential(meters)))
        if villages == 2:
            m1 = 4 * ne * mig[0]  # 4Nm
            return "{}".format(m1)  # mig_matrix is symmetrical and island
        elif villages == 3:
            m1 = 4 * ne * mig[0]
            m2 = 4 * ne * mig[1]
            m3 = 4 * ne * mig[2]
            return "{} {} {} {} {} {} {} {} {}".format(
                0, m1, m2, m1, 0, m3, m2, m3, 0)  # mig_matrix is symmetrical
        elif villages == 4:
            m1 = 4 * ne * mig[0]
            m2 = 4 * ne * mig[1]
            m3 = 4 * ne * mig[2]
            m4 = 4 * ne * mig[3]
            m5 = 4 * ne * mig[4]
            m6 = 4 * ne * mig[5]
            return "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(
                0, m1, m2, m3, m1, 0, m4, m5, m2, m4, 0, m6, m3, m5, m6, 0)

def parse_coalsims_fx(msout, ploidy):
    """Parse output from ms or scrm.
    Parameters
    ---------
    msout: iterator,file
         stdout from scrm
    ploidy: int
         ploidy of the locus
         
    Returns
    ------
    gt_array
         columns to be added to dfAdult
    """
    if ploidy == 1:
         gt_array = []
         for line in iter(msout.stdout.readline, ''):
              if line.startswith("positions"):
                   ##collisions can result here when theta is high
                   pos = np.round(np.array(line.strip().split()[1:], dtype=np.float64))
                   prev = 0
                   for idx, item in enumerate(pos, start=0):
                        while prev >= item:
                             item += 1
                        pos[idx] = item
                        prev = pos[idx]  
                   for line in iter(msout.stdout.readline, ''):
                        hap = np.array(list(line.strip()), dtype=int)
                        gt = hap * pos
                        gt_array.append(gt[gt != 0])
         return gt_array, pos
    elif ploidy == 2:
         gt_array = []
         gt_array2 = []
         for line in iter(msout.stdout.readline, ''):
              if line.startswith("positions"):
                   ##collisions can result here when theta is high
                   pos = np.round(np.array(line.strip().split()[1:], dtype=np.float64))
                   prev = 0
                   for idx, item in enumerate(pos, start=0):
                        while prev >= item:
                             item += 1
                        pos[idx] = item
                        prev = pos[idx]
                   for line in iter(msout.stdout.readline, ''):
                        hap = np.array(list(line.strip()), dtype=int)
                        gt = hap * pos
                        gt_array.append(gt[gt != 0])                     
                        hap2_temp = next(iter(msout.stdout.readline, ''))
                        hap2 = np.array(list(hap2_temp.strip()), dtype=int)
                        gt2 = hap2 * pos
                        gt_array2.append(gt2[gt2 != 0])
         return gt_array, gt_array2, pos         

                   
def coalsims_fx(worm_popsize, villages, initial_migration, initial_distance_m,
        theta, basepairs, mutation, recombination, time2Ancestral, thetaRegional,
        time_join):
    '''External call to ms (Hudson 2000) or scrm.

    Calls the function migration_matrix() if using more than 1 village.
    This function will then read in the stdout from ms/scrm. The
    location of segsites are assigned a random number then scaled
    by the length of the desired sequence in basepairs.

    Parameters
    ----------
    worm_popsize: (list, int)
        how many worms to simulate
    villages: int
        number of villages/metapopulations      
    theta: (list,float)
        list of theta values for each locus. With more than 1 village
        it is list of lists [[locus1_meta1,locus1_meta2],
        [locus2_meta1,locus2_meta2]]
    basepairs:(list,int)
        list of basepairs (lengths) for each locus
    mutation: (list,float)
        mutation rate for each locus as probability per base per generation
    recombination:(list,float)
        recombination rate for each locus as probability per base 
        per generation
    thetaRegional: float
        theta for regional pops to the ratio of N0 e.g.,23 times larger
    time2Ancestral: int
        time in generations to the ancestral population
    time_join: int
        time in generations for joining/splitting villages

    Returns
    -------
    gt : array
         genotype data for first haplotype
    gt2: array
         genotype data for second haplotype 
    pos: array      
         array of mutations positions
    '''
    #theta for first village
    thetaN0 = theta[0]
    #recombination rate for locus
    rho = thetaN0 * (recombination / mutation)    
    rho is 0
    # time to the ancestral population
    tA = basepairs * (time2Ancestral / (thetaN0 / mutation))
    # time to join villages
    tjoin = basepairs * (time_join / (thetaN0 / mutation))
    
    #check ploidy for total haplotypes
    if rho == 0:
         ploidy = 1
    else:
         worm_popsize = [x * 2 for x in worm_popsize]
         ploidy = 2
    #parameters for ms or scrm
    ms_params = {
            'nhaps': sum(worm_popsize),
            'theta': thetaN0,
            'rho': rho,
            'basepair': basepairs - 1,
            'exp_growth': (-1 / tA) * math.log(thetaRegional),
            'time_growth': tA,
            'sig_digits': 12,
        }
     
            #order matters for the command call
    if villages == 1:
        scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                 "-G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                 " {sig_digits} ")
        mscmd = scrm_base.format(**ms_params)
    else:  # ms setup for >1 villages
        num_subpops = len(worm_popsize)  # -I num_pops        
        sub_pop = " ".join(map(str, worm_popsize))  # -I X i j ...
        mm = coalsims_migmat_fx(
           villages,
           initial_migration,
           initial_distance_m,
           thetaN0,
           basepairs,
           mutation)
        if villages == 2:
                 scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                     "{sub_pop} {present_pop} {join} -G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                     " {sig_digits} ")
                 ms_params['present_pop'] = '-n 1 1 -n 2 {}'.format(float(theta[1])/thetaN0)
                 ms_params['sub_pop'] = '-I {} {} {}'.format(num_subpops, sub_pop, mm)
                 ms_params['join'] = '-ej {} 1 2'.format(tjoin)
                 mscmd = scrm_base.format(**ms_params)
        else:     
                 scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                     "{sub_pop} {present_pop} {ma} {join} -G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                     " {sig_digits} ")
                 ms_params['sub_pop'] = '-I {} {}'.format(num_subpops, sub_pop)
                 ms_params['ma'] = '-ma {} '.format(mm)
                 joinstr = ''
                 subpopstr = ''
                 
                 for village_ix in range(villages):
                    present_pop = float(theta[village_ix])/thetaN0
                    subpopstr += '-n {} {} '.format(village_ix + 1, present_pop)
                    if village_ix != villages -1:
                         joinstr += '-ej {0} {1} {2} '.format(
                             tjoin,
                             village_ix + 1,
                             village_ix + 2)
                 ms_params['present_pop'] = subpopstr
                 ms_params['join'] = joinstr
                 mscmd = scrm_base.format(**ms_params)     
    print(mscmd)
    msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)          
    if ploidy == 1:
         gt, mutations = parse_coalsims_fx(msout, ploidy)    
         return gt, mutations
    elif ploidy == 2:
         gt, gt2, mutations = parse_coalsims_fx(msout, ploidy)
         return gt, gt2, mutations
         
def sel_fx(locus, positions):
     '''initializes the distribution of fitness effects for each mutation
     
     Parameters
     ---------
     locus : int
          number of loci
     positions : list
          list from coalsims_fx, the the line "positions" in the scrm/ms output
     Returns
     ------
     dfSel
     '''
     selF =[]
     selS =[]
     positions = [item for sublist in positions for item in sublist]
     #below functs assume the list is nested 
     if isinstance(positions[0], list):
          numpos = [len(i) for i in positions]
          for loc in range(locus):
               for pos in positions[loc]:
                    if random.choice("SF") is "F":
                         #shape = 4, mean = 1, scale = mean/shape
                         #here mean is mean_fitness, wildtype is assumed to be 1
                         selF.append(np.random.gamma(4, scale=0.25))
                         selS.append(1)
                    else:
                         selS.append(np.random.gamma(4, scale=0.25))
                         selF.append(1)
         
          dfSel = pd.DataFrame({
                                'locus' : np.repeat(range(1, locus), numpos),
                                'position' : [item for sub in positions for item in sub],
                                'selF' : selF,
                                'selS' : selS,
                                'freqInt' : np.zeros(sum(numpos))})     
          dfSel = dfSel.loc[:, ['locus', 'position', 'selF',
                 'selS', 'freqInit']]
     else: #list is not nested
          numpos = len(positions)        
          for pos in positions:
               if random.choice("SF") is "F":
                    #shape = 4, mean = 1, scale = mean/shape
                    #here mean is mean_fitness, wildtype is assumed to be 1
                    selF.append(np.random.gamma(4, scale=0.25))
                    selS.append(1)
               else:
                    selS.append(np.random.gamma(4, scale=0.25))
                    selF.append(1)
         
          dfSel = pd.DataFrame({
                                'locus' : np.repeat(range(1, locus), numpos),
                                'position' : [item for item in positions],
                                'selF' : selF,
                                'selS' : selS,
                                'freqInt' : np.zeros(numpos)})     
          dfSel = dfSel.loc[:, ['locus', 'position', 'selF',
                 'selS', 'freqInit']]  
     return dfSel
     
def fit_fx(locus, dfAdult, dfSel):
     ''' calculates mean fitness for each individual by summing fitness effects
     from dfSel for each position across all loci
     
     Parameters
     ---------
     dfAdult : df
          data frame of adult worms containing genotype information
     dfSel : df
          data fram of fitness benefit for each allele
     
     Returns
     ------
     fitF : array
          array filling selF column, influencing fecundity
     fitS : array      
          array filling selS column, influencing survival
     freqinit : array
          initial freq of mutations for dfSel
     '''
     allele = []
    
     ##freq of positions in dfSel
     for loc in range(1, locus):
          for row in range(len(dfAdult)):
               #flatten list and extend
               allele.extend(dfAdult.iloc[row]["locus_" + str(loc) + "_h1"])
               allele.extend(dfAdult.iloc[row]["locus_" + str(loc) + "_h2"])
     freq = np.unique(allele, return_counts=True)[1] / (len(dfAdult) * 2.0) 
     
     ##fitness of individual in dfAdult from values in dfSel               
     fitS_ind = []
     fitF_ind = []
     fitS = []
     fitF = []
     for row in range(len(dfAdult)):
          for loc in range(1, locus):                         
               fitS_ind.extend(dfSel.loc[dfSel["position"].isin(dfAdult.iloc[row]
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selS'][dfSel["locus"] == loc])
               fitF_ind.extend(dfSel.loc[dfSel["position"].isin(dfAdult.iloc[row]
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selF'][dfSel["locus"] == loc])
          fitS.append(round(np.mean(fitS_ind), 5))
          fitF.append(round(np.mean(fitF_ind), 5))  
    
     return fitS, fitF, freq         
               
def wormdf_fx(villages, infhost, muWormBurden, sizeWormBurden, locus,
              initial_migration, initial_distance_m, theta, basepairs, mutation, 
              recombination, time2Ancestral, thetaRegional, time_join, selection):

     #create parasite burden per host
     popinit = []
     for mu, size, numworms in zip(muWormBurden, sizeWormBurden, infhost):
         #the number of worms per infected host
         wb_burden = np.random.negative_binomial(mu, mu/float(mu+size), numworms) 
         #[[10,12,1,15 ... infhostvill1],[19,1,5,4, ... infhostvill2]]
         popinit.append(np.array(wb_burden).tolist())
        
     #total worms per villages   
     wormpopsize = [sum(i) for i in popinit]
     
     host_idx = []
     for vill in range(villages):
          for host in range(len(popinit[vill])):
               host_idx.extend(["v" + str(vill) + "h" + str(host + 1)] * popinit[vill][host])
       
     dfAdult = pd.DataFrame({
                      'village' : np.repeat(range(villages), wormpopsize),
                      'hostidx' : host_idx,
                      'age' : 1,
                      'sex' : [random.choice("MF") for i in range(sum(wormpopsize))],
                      'R0net' : np.random.random(sum(wormpopsize)),
                      'fec' : 0
                      })
     dfAdult = dfAdult.loc[:, ['village', 'hostidx', 'age',
            'sex', 'R0net', 'fec']]
     #add genetic data
     posSel = []
     for loc in range(locus):
          if recombination[loc] == 0:
               gt_array, mutations = coalsims_fx(wormpopsize, villages, initial_migration, 
                                    initial_distance_m, theta[loc], basepairs[loc], mutation[loc], 
                                    recombination[loc], time2Ancestral, thetaRegional, time_join)
               #gt_array, mutations = coalsims_fx([10],1,.0001,[0],[5],13000,7.6E-8,0,1800,23,240)
               #gt_array, mutations = coalsims_fx([10,10],2,.0001,[1000],[5,5],13000,7.6E-8,0,1800,23,240)
               dfAdult["locus_" + str(loc)] = gt_array          
          elif recombination[loc] > 0:
               gt_array, gt_array2, mutations = coalsims_fx(wormpopsize, villages, initial_migration, 
                          initial_distance_m, theta[loc], basepairs[loc], mutation[loc], 
                          recombination[loc], time2Ancestral, thetaRegional, time_join)
               #gt_array, gt_array2, mutations = coalsims_fx([10],1,.0001,[0],[5],13000,7.6E-8,2.9E-9,1800,23,240)
               #gt_array, gt_array2, mutations = coalsims_fx([10,10],2,.0001,[1000],[5,5],13000,7.6E-8,2.9E-9,1800,23,240) 
               dfAdult["locus_" + str(loc) + "_h1"] = gt_array
               dfAdult["locus_" + str(loc) + "_h2"] = gt_array2
               posSel.append(mutations)                      

     #create dfSel
     if selection:
          dfSel = sel_fx(locus, posSel)
          fitS, fitF, freq = fit_fx(locus, dfAdult, dfSel)     
          dfSel["freqInit"] = freq                   
          dfAdult["fitF"] = fitF
          dfAdult["fitS"] = fitS          
          return dfAdult, dfSel
     else:
          return dfAdult

def wbsims_init(villages, hostpopsize, prevalence, muTrans, sizeTrans, muWormBurden, 
                sizeWormBurden, locus, initial_migration, initial_distance_m, theta,
                basepairs, mutation, recombination, time2Ancestral, thetaRegional,
                time_join, selection):
     '''main function call for simulations
     Parameters
     ---------
     villages : int
          number of villages to simulate
     villpopulation : int, list
          maximum population size of each village
     prevalence : float, list
          prevalence of disease in each village
     
          
     Returns
     ------
     dfAdult : df
          initialize df for adult worms
     dfHost : df
          initialize df for host
     dfSel :df 
          initialized df for selection
     '''
     hostpopsize = np.array(hostpopsize)
     prevalence = np.array(prevalence)
     infhost = np.round(hostpopsize * prevalence).astype(np.int64)
     #construct dfHost
     dfHost = host_fx(villages, infhost, muTrans, sizeTrans)
     # dfHost = host_fx(2, [10,10], 100, 1)
     #construct dfAdult, dfSel
     if selection:
          dfAdult, dfSel= wormdf_fx(villages, infhost, muWormBurden, sizeWormBurden, 
                               locus, initial_migration, initial_distance_m, theta, 
                               basepairs, mutation, recombination, time2Ancestral, thetaRegional,
                               time_join, selection)
     else:
          dfAdult = wormdf_fx(villages, infhost, muWormBurden, sizeWormBurden, 
                     locus, initial_migration, initial_distance_m, theta, 
                     basepairs, mutation, recombination, time2Ancestral, thetaRegional,
                     time_join, selection)
     #dfAdult = wormdf_fx(2, [10,10], [5,5], [50,50], 1, 0.0001, [1000], [[5,5],[]], [13000], 
     #                            [7.6E-8], [0], 1800, 23, 240, False)
     
     #dfAdult, dfSel = wormdf_fx(2, [10,10], 5, 50, 2, 0.0001, [1000], [[5,5],[10,10]], [13000,200000], 
     #                            [7.6E-8, 2.9E-9], [0, 2.9E-9], 1800, 23, 240, True)
     dfJuv = pd.DataFrame({})
     dfMF = pd.DataFrame({})
     if selection:
          return dfAdult, dfHost, dfSel, dfJuv, dfMF  
     else:
          return dfAdult, dfHost, dfJuv, dfMF
     
def hostmigration_fx(dfHost):
     ''' allows host to move between villages
     Parameters
     ---------
     dfHost : df
          df of hosts
     
     Returns
     --------
     dfHost : df
     ''' 
     migrant = np.random.randint(0, len(dfHost))
     if dfHost.loc[migrant, "village"] < max(dfHost.village):
          dfHost.loc[migrant, "village"] += 1
     else:
          dfHost.loc[migrant, "village"] -= 1                         

     return dfHost

def vectorbite_fx(bitespperson, 
                  hours2bite, 
                  hostpopsize,
                  prev_t,
                  densitydep, 
                  avgMF, 
                  bednets, 
                  bnstart, 
                  bnstop, 
                  bncoverage, 
                  month):
    '''counts the number of successful infectious mosquito bites
    Parameters
    --------
    bitespperson : int
        rate of bites per person per unit time, here hours
    hours2bite : int
        number of hours mosq bite per night/day
    hostpopsize : int
        total population of village 
    prev_t : float
        percent of people infected 
    densityDep : Boolean
        use the density dependece function for developing L3
    avgMF : float
         average number of MF per host per village
    bednets : Boolean, list
        use values from bednets     
    bnstart : int, list
        month to start bednets 
    bnstop : int, list
        mont to stop bednets  
    bncoverage : int, list
        percent of population using bednets 
    month : int
        current time, month 
    Returns
    ------
    L3trans : int
         all new juv are age class 0
    '''
    print "vectorbite"
    if bednets:
        if month > bnstart and month < bnstop:
             totalbites = ((1 - bncoverage) * bitespperson * hours2bite * 30 
                                   * hostpopsize)
        else:
             totalbites = (bitespperson * hours2bite * 30 
                                   * hostpopsize)
    else:
        totalbites = (bitespperson * hours2bite * 30 
                              * hostpopsize)
    # 0.37 is prob of bite on infected host picking up MF    
    infbites = np.random.binomial(totalbites, (prev_t * 0.37))
    if densitydep: #values for anopheles from CITE
       #number of MF in 20ul of blood
       #235ml is 5% of total host blood 
       #there are 50 units of 20ul in 1ml
       mfBlood = avgMF / 50.0 
       # 0.414 is proportion of  L3 that leave mosquito per bite 
       # 0.32 proportion of L3 that make it into the host
       L3trans = round(infbites * (4.395 * (1 - math.exp( -(0.055 * (mfBlood)) 
           / 4.395)) ** 2) * (0.414 * 0.32)) 
    else:
       # 0.414 is proportion of  L3 that leave mosquito per bite 
       # 0.32 proportion of L3 that make it into the host
       L3trans = np.random.binomial(infbites, (0.414 * 0.32))
    print int(L3trans)
    return int(L3trans)

def new_infection_fx(dispersal, transMF, disthost, dfHost, deathdict):
    '''a transmission event infecting a naive host
    Parameters
    ----------
    mindist : int
         minimum distance
    dispersal: float
         dispersal distance as 2*sigma
    dfHost: df
        host dataframe 
    transMF
        row of transmitted MF
    Returns
    ------
    dfHost: df
         updated with new host information
    newhostidx: str
          new host index     
    '''
    print "newfx"
    #how close
    mindist = 5
    #how far
    maxdist = dispersal
    #assign new position
    maxc = int( math.sqrt(( maxdist ** 2 / 2 )))
    #makes the below statement true
    newpts = dfHost[dfHost.hostidx == transMF.hostidx].coordinates.values[0]

    while next((True for elem in dfHost['coordinates'].values if np.array_equal(elem, newpts)),False) is True:
         #shift of new points
         print "loop"
         x_new = random.randint( -maxc, maxc )
         y_new = random.randint( -maxc, maxc )
         #test that new shifts are min distance
         if math.sqrt(x_new ** 2 + y_new ** 2) > mindist:
              newhostptX = dfHost.coordinates[dfHost.hostidx == transMF.hostidx].values[0][0] + x_new
              newhostptY = dfHost.coordinates[dfHost.hostidx == transMF.hostidx].values[0][1] + y_new
              newpts = np.array([newhostptX, newhostptY])
              print newpts
    print "outloop"          
    #copy village
    vill = transMF.village
    #new host index
    old_hostidx = dfHost[dfHost.village == transMF.village].hostidx.iloc[-1]
    new_hostidx = old_hostidx[:old_hostidx.rfind('h')] + 'h' + str(int(old_hostidx.split('h')[1]) + 1)
    #radnom sex
    sex = random.choice("01")
    #age and agedeath function from wbsims_initialize
    age, agedeath = agehost_fx(sex, deathdict)
    #add to dfHost at bottom
    dfHost.loc[len(dfHost) + 1] = [vill, new_hostidx, sex, age, agedeath, newpts, 0]
    
    return dfHost, new_hostidx

def transmission_fx(villages, hostpopsize, sigma, bitesPperson, hours2bite, densityDep, bednets,
                    bnstart, bnstop, bncoverage, month, dfMF, dfJuv, dfHost, deathdict):
     #transmission_fx(2, [200, 200], 50, [20, 20], [8, 8], [True, True], [False, False], 
                      #[12, 12], [36, 36], 1, dfMF, dfJuv, dfHost, deathdict) 
    '''transmission events resolved as either reinfection or new infection
    Parameters
    --------
    villages : int
        number of villages
    hostpopsize : int, list
        total population size of village      
    sigma : float
        dispersal mean
    bitespperson : int, list
        rate of bites per person per unit time, here hours
    hours2bite : int, list
        number of hours mosq bite per night/day
    densityDep : Boolean
        use the density dependece function for developing L3
    bednets : Boolean, list
        use values from bednets     
    bnstart : int, list
        month to start bednets 
    bnstop : int, list
        mont to stop bednets  
    bncoverage : int, list
        percent of population using bednets 
    month : int
        current time, month 
    dfHost: df
         chooses the donating and accepting hosts for transmission
    dfMF: df
         donating MF genotype to transmit
    dfJuv: df
         donated MF move to Juvenille age class
    Returns
    ------
    dfJuv : df
         all new juv are age class 0
    dfHost : df
         in the case of newinfection, new host
    '''
    print "transmission"
    dispersal = 2 * sigma
    for vill in range(villages):
         prev_t = len(dfHost[dfHost.village == vill]) / float(hostpopsize[vill])
         infhost = len(dfHost[dfHost.village == vill])
         avgMF = len(dfMF[dfMF.village == vill])/float(infhost)
         avgMF = 200
         L3trans = vectorbite_fx(bitesPperson[vill], hours2bite[vill], hostpopsize[vill], 
                                 prev_t, densityDep[vill], avgMF, bednets[vill], 
                                 bnstart[vill], bnstop[vill], bncoverage[vill], month)    
         if L3trans > len(dfMF[dfMF.village == vill]):
              transMF = dfMF[dfMF.village == vill]     
         else:
              transMF = dfMF[dfMF.village == vill].sample(L3trans)
         distMat = pairwise_distances(np.vstack(dfHost[dfHost.village == vill].coordinates))
###alternative way only returns points whose distance is <= dispersal
         #meaning these are the only hosts that can transmit to each other and 
         #disMat is a dict with keys as tuples
         #from scipy.spatial import KDTree
         #tree = KDTree(np.vstack(dfHost[dfHost.village == vill].coordinates))
         #distMat = tree.sparse_distance_matrix(np.vstack(dfHost[dfHost.village == vill].coordinates), dispersal)
         #rehost = [i for i, x in enumerate(distMat.keys()) if x[0] == index]                                    
               
         for index, row in transMF.iterrows():
              dfdistHost = dfHost[dfHost.village == vill]
              #new infection
              disthost = np.where((dfdistHost["hostidx"] == row.hostidx))[0]
              if len(dfHost[dfHost.village == vill]) < hostpopsize[vill]:
                   prob_newinfection = 1.0/(len(distMat[disthost] <= dispersal) + 1)
              else: #everyone is already infected
                   prob_newinfection = 0
              
              print prob_newinfection, row, index 
              if np.random.random() < prob_newinfection:     
                   print "new loop"
                   #new host
                   dfHost, newidx = new_infection_fx(dispersal, row, disthost, dfHost, deathdict)
                   row.hostidx = newidx
                   row.age = 0
                   dfJuv = dfJuv.append(row)
                   #need to update distMat to include new host
                   distMat = pairwise_distances(np.vstack(dfHost[dfHost.village == vill].coordinates)) 
              else: #reinfection
                   print "reinfect loop"
                   rehost = dfHost.iloc[random.choice(np.where((distMat[disthost] <= dispersal)[0])[0])]
                   row.hostidx = rehost.hostidx
                   row.age = 0
                   dfJuv = dfJuv.append(row)                          
    return dfHost, dfJuv

def mutation_fx(locus, dfAdult_mf, mutation_rate, recombination_rate, basepairs):
    '''calculates number of mutations, changes haplotypes
    
    Parameters
    ---------
    locus: int
         number of loci
    dfAdult_mf : pandas dataframe
          dataframe containing new larval parasites
    mutation_rates : float, list
          mutation rates for each locus 
    basepairs : int, list
          length of each locus in basepairs       
    
    Returns
    ------
    dfAdult_mf : pandas df
         df of larval parasites with mutated positions
    mutations : int, list
         list of positions of new mutations
    '''
    positions = []
    muts_counter = []      
    for loc in range(locus):
         if recombination_rate[loc] == 0:
              num_muts = np.random.binomial(len(dfAdult_mf), basepairs[loc] * mutation_rate[loc])
              #print num_muts
              muts_counter.append(num_muts)
              muts = 0     
              if num_muts != 0:
                   while muts < num_muts:
                        randmf = np.random.randint(0,len(dfAdult_mf))
                        newsite = np.random.randint(1,basepairs[loc])
                        positions.append(newsite)
                        newhap = np.append(dfAdult_mf.iloc[randmf]["locus_" + str(loc)], newsite)
                        dfAdult_mf.set_value(randmf, "locus_" + str(loc), newhap) 
                        muts += 1
          
         else:
              num_muts = np.random.binomial(2 * len(dfAdult_mf), basepairs[loc] * mutation_rate[loc])
              #print num_muts
              muts_counter.append(num_muts)
              muts = 0     
              if num_muts != 0:
                   while muts < num_muts:
                        randmf = np.random.randint(0,len(dfAdult_mf))
                        newsite = np.random.randint(1,basepairs[loc])
                        positions.append(newsite)
                        randhap = random.choice("12")
                        newhap = np.append(dfAdult_mf.iloc[randmf]["locus_" + str(loc) + "_h" + randhap], newsite)
                        dfAdult_mf.set_value(randmf, "locus_" + str(loc) + "_h" + randhap, newhap.sort()) 
                        muts += 1
    #print positions
    #print muts_counter                    
    dfMuts = pd.DataFrame({
                           "locus" : np.repeat(range(locus), muts_counter),
                           "position" : positions,
                           "selF" :  np.zeros(len(positions)),
                           "selS" :  np.zeros(len(positions)),
                           "freqInit" :  np.zeros(len(positions))
                           })
    return dfAdult_mf, dfMuts
 
def recombination_fx(locus, dfAdult, dfAdult_mf, recombination_rate, basepairs):
    """calculate number of recombination events and rearranges haplotypes

    Parameters
    ---------
    locus: int
         number of loci
    dfAdult_mf : pandas dataframe
          dataframe containing new larval parasites
    dfAdult : pd df
          dataframe containing reproducing adults      
    recombination_rate : float, list
          recombination rate for each locus 
    basepairs : int, list
          length of each locus in basepairs 

    Returns
    -------
    dfAdult_mf : pd df
    
    """
    dfAdult_mf = pd.DataFrame({})

    for index, row in dfAdult[dfAdult.sex == "F"].iterrows(): #for each female
         try:
              male = dfAdult.loc[(dfAdult["sex"] == "M") & (dfAdult["hostidx"] == row.hostidx)].sample(1)
         except ValueError:
              print "no males, no sex"
              break
         mf = 0
         while mf < dfAdult.loc[index, "fec"]:
              for loc in range(locus):
                   if recombination_rate[loc] == 0:
                        pass
                   else:                                         
                        num_recomb = np.random.poisson(recombination_rate[loc] * basepairs[loc] * 2)
                        print num_recomb
                        if num_recomb == 0:
                             row["locus_" + str(loc) + "_h1"] = row["locus_" + str(loc) + "_h" + random.choice("12")]     
                             #male contribution
                             row["locus_" + str(loc) + "_h2"] = male["locus_" + str(loc) + "_h" + random.choice("12")]                                                    
                        else:
                             r = 0
                             #randomly choose male or female
                             sex_xing = random.choice("MF")
                             #while loop to account for multiple recombination events
                             h1m = np.concatenate(male["locus_" + str(loc) + "_h1"].values)
                             h2m = np.concatenate(male["locus_" + str(loc) + "_h2"].values)
                             h1f = np.concatenate(row["locus_" + str(loc) + "_h1"].values)
                             h2f = np.concatenate(row["locus_" + str(loc) + "_h2"].values)
                             if sex_xing is "M":
                                  while r < num_recomb:
                                       crossover_pos = random.randint(0, basepairs[loc])
                                       try:
                                            hap1_co = np.where(h1m > crossover_pos)[0][-1]
                                            hap2_co = np.where(h2m > crossover_pos)[0][-1]
                                       except IndexError:
                                            break
                                       h1m_new = h1m[0:hap1_co + 1] + h2m[hap2_co:] 
                                       h2m_new = h2m[0:hap2_co + 1] + h1m[hap1_co:]
                                       h1m = h1m_new
                                       h2m = h2m_new
                                       r += 1
                             elif sex_xing is "F":
                                  while r < num_recomb:
                                       crossover_pos = random.randint(0, basepairs[loc])
                                       try:
                                            hap1_co = np.where(h1f> crossover_pos)[0][-1]
                                            hap2_co = np.where(h2f > crossover_pos)[0][-1]
                                       except IndexError:
                                            break
                                       h1f_new = h1f[0:hap1_co + 1] + h2f[hap2_co:] 
                                       h2f_new = h2f[0:hap2_co + 1] + h1f[hap1_co:]
                                       h1f = h1f_new
                                       h2f = h2f_new
                                       r += 1
                             row["locus_" + str(loc) + "_h1"] = random.choice([h1f, h2f])     
                             row["locus_" + str(loc) + "_h2"] = random.choice([h1m, h2m])      
              dfAdult_mf = dfAdult_mf.append([row], ignore_index=True)
              mf += 1
     
    return dfAdult_mf

def selection_fx(dfAdult_mf, dfMuts, dfSel, locus):
     '''recalculates DFE for new mutations and phenotype for new mf
     Parameters
     ---------
     dfSel : df
          updates this dataframe
     dfAdult_mf : df
          adds phenotype info
     dfMuts : df
          gives positions of new mutations for addition to dfSel
     locus : int
          number of loci     
     Returns
     ------
     dfSel : df
          now updates with new mutations
     dfAdult_mf : df
          updated with phenotype     
     '''
     for loc in range(1, locus):
          for index, row in dfMuts.iterrows():
               if random.choice("SF") is "F":
                    #shape = 4, mean = 1, scale = mean/shape
                    #here mean is mean_fitness, wildtype is assumed to be 1
                    row.selF.set_value(np.random.gamma(4, scale=0.25))
                    row.selS.set_value(1)
               else:
                    row.selS.set_value(np.random.gamma(4, scale=0.25))
                    row.selF.set_value(1)        
     dfSel = dfSel.append(dfMuts)
     dfAdult_mf = fitness_fx(locus, dfAdult_mf, dfSel)               
     return dfAdult_mf, dfSel
     
def fitness_fx(locus, dfAdult_mf, dfSel):
     ''' calculates mean fitness for each individual by summing fitness effects
     from dfSel for each position across all loci
     
     Parameters
     ---------
     dfAdult_mf : df
          data frame of adult worms containing genotype information
     dfSel : df
          data fram of fitness benefit for each allele
     
     Returns
     ------
     dfAdult_mf : df
          updated df for mf
     '''
     ##fitness of individual in dfAdult from values in dfSel               
     for index, row in dfAdult_mf.iterrows():
          fitS_ind = []
          fitF_ind = []
          for loc in range(1,locus):   
               fitS_ind.extend(dfSel.loc[dfSel["position"].isin(row
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selS'][dfSel["locus"] == loc])
               fitF_ind.extend(dfSel.loc[dfSel["position"].isin(row
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selF'][dfSel["locus"] == loc])
          row.fitS.set_value(round(np.mean(fitS_ind), 5))
          row.fitF.set_value(round(np.mean(fitF_ind), 5))  
    
     return dfAdult_mf

def fecunditybase_fx(fecund, dfAdult, locus, mutation_rate, recombination_rate, basepairs, selection):
    #(20, dfAdult, 2, [7.6E-8, 2.9E-9], [0, 2.9E-9], [13000, 200000])
     '''base fecundity function, simpliest scenario
    conditions: mda=False, selection=F
    
    Parameters
    ---------
    fecund: int
         rate of the poisson distribution for offspring number
    dfAdult: pandas dataframe
          dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites 
          
    Returns
    ------
    dfAdult_mf : df
         deep copy of adult genotypes
    
    TO DO: track the rate of decline for fecundity, should be linear
     
     '''
     #all locations where age is less than 6
     dfAdult.loc[dfAdult.age < 6, "fec"] = np.random.poisson(fecund, 
               len(dfAdult[dfAdult.age < 6]))         
     #linear function defining decline in fecundity with age
     m = float(0 - fecund) / (21 - 6)
     b = 0 - m * 21
     #assign fecundity value based on age function     
     dfAdult.loc[dfAdult.age >= 6, "fec"] = np.random.poisson(m 
               * dfAdult.loc[dfAdult.age >= 6,"age"] + b)

     #sex, recombination, mutation
     dfAdult_mf = recombination_fx(locus, dfAdult, recombination_rate, basepairs)
     dfAdult_mf, dfMuts = mutation_fx(locus, dfAdult_mf, mutation_rate, recombination_rate, 
                                      basepairs, selection)
    
     if selection:
         dfAdult_mf, dfSel = selection_fx(dfAdult_mf, dfMuts, dfSel)
    
     return dfAdult_mf, dfSel if selection is True else dfAdult_mf

def survivalbase_fx(month, surv_Juv, shapeMF, scaleMF, shapeAdult, scaleAdult,
                    dfMF, dfAdult, dfJuv, dfHost, fecund, locus, recombination_rate,
                    mutation_rate, basepairs, selection):
   #(1, 0.866, 3.3, 10, 3.8, 8, dfMF, dfAdult, dfJuv, dfHost)  
                        
    '''base survival function
    Parameters
    ---------
    month: int
         time in months of the simulation
    surv_Juv: float
         survival prob of juvenille life stage
    shapeMF: float
         shape parameter for weibull distribution of MF
    scaleMF: int
         scale parameter for weibull distribution of MF
    shapeAdult: float
         shape parameter for weibull distribution of Adults
    scaleAdult: int
         scale parameter for weibull distribution of MF

    Returns
    -------
    dfMF
    dfAdult
    dfJuv
    dfHost
    
    '''
    #adult worms and hosts are only evaluated per year
    if month%12 == 0:
        #Adult survival is based on weibull cdf
        surv_adultrand = np.random.random(len(dfAdult))    
        surv_adultfxage = weibull_min.cdf(dfAdult.age, shapeAdult,loc=0,scale=scaleAdult)
        surviveAdult = np.where(surv_adultrand <= (1 - surv_adultfxage))
        dfAdult = dfAdult.iloc[surviveAdult]
        dfAdult.age = dfAdult.age + 1 #2 - 21
        
        ##host survival is from act table
        dfHost = dfHost[dfHost.age < dfHost.agedeath]
        #remove all worms with dead host.hostidx from all dataframes
        dfAdult = dfAdult.loc[dfAdult["hostidx"].isin(dfHost.hostidx)]
        dfJuv = dfJuv.loc[dfJuv["hostidx"].isin(dfHost.hostidx)]
        dfMF = dfMF.loc[dfMF["hostidx"].isin(dfHost.hostidx)]
        #add 1 year to all ages of hosts
        dfHost.age = dfHost.age + 1 
        
    ##Juv is exponential 0.866; surv_Juv
    #dont include age 0 which just moved from transmission fx
    surv_juvrand = np.random.random(len(dfJuv.age > 0))
    surviveJuv = np.where(surv_juvrand <= surv_Juv)
    dfJuv = dfJuv.iloc[surviveJuv]
    dfJuv.age = dfJuv.age + 1 # 1 - 13

    ##MF is weibull cdf
    surv_mfrand = np.random.random(len(dfMF))    
    surv_mffxage = weibull_min.cdf(dfMF.age,shapeMF,loc=0,scale=scaleMF)
    surviveMF = np.where(surv_mfrand <= (1 - surv_mffxage))
    dfMF = dfMF.iloc[surviveMF]
    dfMF.age = dfMF.age + 1 #2 - 12
    dfMF = dfMF[dfMF.age < 13] #hard cutoff at 12 months

    ##move Juv age 13 to adult age 1
    #dfJuv_new = pd.DataFrame({})
    dfJuv_new = dfJuv[dfJuv.age > 12]
    #reset age to adult   
    dfJuv_new.age = 1
    #increase R0net for next gen
    dfJuv_new.R0net = dfJuv_new.R0net + 1
    #append to adults
    dfAdult = dfAdult.append(dfJuv_new, ignore_index=True)
    #remove Juv age 13 from dfJuv
    dfJuv = dfJuv[dfJuv.age <= 12]
     
    ##call to fecundity fx to deepcopy adult to dfMF age 1
    #fecundity calls mutation/recombination
    dfAdult_mf = fecunditybase_fx(fecund, dfAdult, locus, mutation_rate, recombination_rate, 
                                  basepairs, selection)
    dfAdult_mf.age = 1
    dfAdult_mf.fec = 0
    dfAdult_mf.sex = [random.choice("MF") for i in range(len(dfAdult_mf))]
    dfMF = dfMF.append(dfAdult_mf, ignore_index=True)     

    return dfAdult, dfHost, dfJuv, dfMF if month%12 == 0 else dfJuv, dfMF
 
def wb_sims(numberGens):
    '''main function for simulations
    Parameters
    ---------
    numberGens : int
         how many months to run the simulation
    burning : int
         burning before recording data
    config_file : file
         file with options to be read by configparser()
    Returns
    -------
    
    '''
    #wbsims_init
    villages = 2
    hostpopsize = [100, 200]
    prevalence = [0.1, 0.3]
    muTrans = 100
    sizeTrans = 1
    muWormBurden = [5, 5]
    sizeWormBurden = [50, 50]
    locus = 2
    initial_migration = 0.0001
    initial_distance_m = [1000]
    theta = [[5, 5], [10, 10]]
    basepairs = [13000, 200000]
    mutation_rate = [7.6E-8, 2.9E-9]
    recombination_rate = [0, 2.9E-9]
    time2Ancestral = 1800
    thetaRegional = 23
    time_join = 240
    selection = False
    #transmission
    sigma = 100
    bitesPperson = [10, 10]
    hours2bite = [8, 8]
    densityDep = [True, True]
    bednets = [False, False]
    bnstart = [0, 0]
    bnstop = [0, 0]
    bncoverage = [0, 0]
    #survival
    surv_Juv = 0.866
    shapeMF = 3.3
    scaleMF = 10
    shapeAdult = 3.8
    scaleAdult = 8
    fecund = 20
    
    #load needed actuarial table from file
    deathdict = {}
    with open("/home/scott/act.tbl",'r') as tbl:
         for line in tbl:
              line = line.strip()
              deathdict["{}".format(line.split()[0])] = list(map(float,
                  line.split()[1:]))
    #set counters
    month = 0
    sim_time = numberGens

    if selection:
         dfAdult, dfHost, dfMF, dfJuv, dfSel = wbsims_init(villages, hostpopsize, prevalence, muTrans, sizeTrans, 
                                  muWormBurden, sizeWormBurden, locus, initial_migration, 
                                  initial_distance_m, theta, basepairs, mutation_rate, 
                                  recombination_rate, time2Ancestral, thetaRegional,
                                  time_join, selection)
         while month <= sim_time:
              dfJuv, dfHost = transmission_fx(villages, hostpopsize, sigma, bitesPperson, 
                                        hours2bite, densityDep, bednets, bnstart,
                                        bnstop, bncoverage, month, dfMF, dfJuv, dfHost, deathdict)
              dfAdult, dfJuv, dfMF, dfHost, dfSel = survivalbase_fx(month, surv_Juv, shapeMF, scaleMF, shapeAdult,
                                                   scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                                                   fecund, locus, mutation_rate, recombination_rate, 
                                                   basepairs, selection, dfSel) 
              month += 1

         
    else:
         dfAdult, dfHost, dfMF, dfJuv = wbsims_init(villages, hostpopsize, prevalence, muTrans, sizeTrans, 
                                  muWormBurden, sizeWormBurden, locus, initial_migration, 
                                  initial_distance_m, theta, basepairs, mutation_rate, 
                                  recombination_rate, time2Ancestral, thetaRegional,
                                  time_join, selection)        
         while month <= sim_time:
             dfJuv, dfHost = transmission_fx(villages, hostpopsize, sigma, bitesPperson, 
                                             hours2bite, densityDep, bednets, bnstart,
                                             bnstop, bncoverage, month, dfMF, dfJuv, dfHost, deathdict)
             dfAdult, dfJuv, dfMF, dfHost = survivalbase_fx(month, surv_Juv, shapeMF, scaleMF, shapeAdult,
                                                        scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                                                        fecund, locus, mutation_rate, recombination_rate, 
                                                        basepairs, selection) 
             month += 1

if __name__ == '__main__':
     wb_sims(24)
 