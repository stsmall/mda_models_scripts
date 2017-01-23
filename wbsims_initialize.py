import math
import subprocess
import numpy as np
import pandas as pd
import random

random.seed(2000)
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
    with open("/home/scott/act.tbl",'r') as tbl:
         for line in tbl:
              line = line.strip()
              deathdict["{}".format(line.split()[0])] = list(map(float,
                  line.split()[1:]))
    
    #when do they die
    death = deathdict[str(age)][int(sex)] + np.random.normal(0,6)
         
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
    
def migration_matrix(villages, initial_migration, initial_distance_m, theta,
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

def parse_ms_output(msout, ploidy):
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
                   for line in iter(msout.stdout.readline, ''):
                        hap = np.array(list(line.strip()), dtype=int)
                        gt = hap * pos
                        gt_array.append(gt[gt != 0])                     
                        hap2_temp = next(iter(msout.stdout.readline, ''))
                        hap2 = np.array(list(hap2_temp.strip()), dtype=int)
                        gt2 = hap2 * pos
                        gt_array2.append(gt2[gt2 != 0])
         return gt_array, gt_array2, pos         

                   
def ms_outcall(worm_popsize, villages, initial_migration, initial_distance_m,
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
        mm = migration_matrix(
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
         gt, mutations = parse_ms_output(msout, ploidy)    
         return gt, mutations
    elif ploidy == 2:
         gt, gt2, mutations = parse_ms_output(msout, ploidy)
         return gt, gt2, mutations
         
def sel_fx(locus, positions):
     '''initializes the distribution of fitness effects for each mutation
     
     Parameters
     ---------
     locus : int
          number of loci
     positions : list
          list from ms_outcall, the the line "positions" in the scrm/ms output
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
     
def fitness_fx(locus, dfAdult, dfSel):
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
               gt_array, mutations = ms_outcall(wormpopsize, villages, initial_migration, 
                                    initial_distance_m, theta[loc], basepairs[loc], mutation[loc], 
                                    recombination[loc], time2Ancestral, thetaRegional, time_join)
               #gt_array, mutations = ms_outcall([10],1,.0001,[0],[5],13000,7.6E-8,0,1800,23,240)
               #gt_array, mutations = ms_outcall([10,10],2,.0001,[1000],[5,5],13000,7.6E-8,0,1800,23,240)
               dfAdult["locus_" + str(loc)] = gt_array          
          elif recombination[loc] > 0:
               gt_array, gt_array2, mutations = ms_outcall(wormpopsize, villages, initial_migration, 
                          initial_distance_m, theta[loc], basepairs[loc], mutation[loc], 
                          recombination[loc], time2Ancestral, thetaRegional, time_join)
               #gt_array, gt_array2, mutations = ms_outcall([10],1,.0001,[0],[5],13000,7.6E-8,2.9E-9,1800,23,240)
               #gt_array, gt_array2, mutations = ms_outcall([10,10],2,.0001,[1000],[5,5],13000,7.6E-8,2.9E-9,1800,23,240) 
               dfAdult["locus_" + str(loc) + "_h1"] = gt_array
               dfAdult["locus_" + str(loc) + "_h2"] = gt_array2
               posSel.append(mutations)                      

     #create dfSel
     if selection:
          dfSel = sel_fx(locus, posSel)
          fitS, fitF, freq = fitness_fx(locus, dfAdult, dfSel)     
          dfSel["freqInit"] = freq                   
          dfAdult["fitF"] = fitF
          dfAdult["fitS"] = fitS          
          return dfAdult, dfSel
     else:
          return dfAdult, posSel

def wbsims_init(villages, villpopulation, prevalence, muTrans, sizeTrans, muWormBurden, 
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
     villpopulation = np.array(villpopulation)
     prevalence = np.array(prevalence)
     infhost = np.round(villpopulation * prevalence).astype(np.int64)
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
     if selection:
          return dfAdult, dfHost, dfSel  
     else:
          return dfAdult, dfHost

if __name__ == '__main__':
     #2 villages with selection
     dfAdult, dfHost, dfSel = wbsims_init(2, [100, 200], [0.1, 0.3], 100, 1, [5, 5], [50, 50], 2, 0.0001, [1000], 
                               [[5, 5], [1, 1]], [13000, 200000], 
                               [7.6E-8, 2.9E-9], [0, 2.9E-9], 1800, 23, 240, True)    
     #2 villages without selection
     dfAdult, dfHost = wbsims_init(2, [100, 200], [0.1, 0.3], 100, 1, [5, 5], [50, 50], 2, 0.0001, [1000], 
                               [[5, 5], [1, 1]], [13000, 200000], 
                               [7.6E-8, 2.9E-9], [0, 2.9E-9], 1800, 23, 240, False) 
