import math
import subprocess
import numpy as np
import pandas as pd
import random

from IPython import embed

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
    with open("act.tbl",'r') as tbl:
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
             host_idx.append("v" + str(vill) + "h" + str(host))
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
    gt_array = []
    for line in iter(msout.stdout.readline, ''):
         if line.startswith("positions"):
              pos = np.round(np.array(line.strip().split(), dtype=np.float64))
              for line in iter(msout.stdout.readline, ''):
                   hap = np.array(list(line.strip()), dtype=int)
                   if ploidy is 1: #haploid
                        gt_array = pos * hap
                   else: #diploid
                        hap2 = np.array(list(next(line).strip()), dtype=int)
                        gt_array = np.vstack(pos * hap, pos * hap2)
    return gt_array, pos 

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
         genotype data 
    '''
    #theta for first village
    thetaN0 = theta[0]
    #recombination rate for locus
    rho = thetaN0 * (recombination / mutation)    
    # time to the ancestral population
    tA = basepairs * (time2Ancestral / (thetaN0 / mutation))
    # time to join villages
    tjoin = basepairs * (time_join / (thetaN0 / mutation))
    
    #check ploidy for total haplotypes
    if rho is 0:
         total_inds = sum(worm_popsize)
         ploidy = 1
    else:
         total_inds = 2 * sum(worm_popsize)
         ploidy = 2
    #parameters for ms or scrm
    ms_params = {
            'nhaps': total_inds,
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
    gt, mutations = parse_ms_output(msout, ploidy)    
          
    return gt, mutations
    
def sel_fx(positions, seloption):
     '''initializes the distribution of fitness effects for each mutation
     
     Parameters
     ---------
     positions: list
          list from ms_outcall, the the line "positions" in the scrm/ms output
     seloption: int
          0 is no selection
          1 is no historic seleciton
          2 is historic selection
     Returns
     ------
     dfSel
     '''
     dfSel = pd.DataFrame({
                           'locus' : [],
                           'position' : [],
                           'selF' : [],
                           'selS' : [],
                           'freqInt' : []})     
     if seloption is 0:
          pass
          return
     elif seloption is 1: #neutral to history, DFE from gamma
          idx = 0
          for loc in range(len(positions)):  #2 loci would be 0,1
               for pos in positions[loc]:
                    dfSel[idx]["locus"] = "locus_" + str(loc)
                    dfSel[idx]["position"] = positions[loc][pos]
                    if random.choice("SF") is "F":
                         dfSel[idx]["selF"] = np.random.gamma()
                         dfSel[idx]["selS"] = 0
                    else:
                         dfSel[idx]["selS"] = np.random.gamma()
                         dfSel[idx]["selF"] = 0                         
                    idx += 1
          return dfSel
     elif seloption is 2: #non-neutral to history, DFE from normal
          idx = 0
          for loc in range(len(positions)):  #2 loci would be 0,1
               for pos in positions[loc]:
                    dfSel[idx]["locus"] = "locus_" + str(loc)
                    dfSel[idx]["position"] = positions[loc][pos]
                    if random.choice("SF") is "F":
                         dfSel[idx]["selF"] = np.random.norm()
                         dfSel[idx]["selS"] = 0
                    else:
                         dfSel[idx]["selS"] = np.random.norm()
                         dfSel[idx]["selF"] = 0                         
                    idx += 1     
          return dfSel
     
def fitness_fx(dfAdult, dfSel):
     ''' calculates fitness for each individual by summing fitness effects
     from dfSel for each position across all loci
     
     Parameters
     ---------
     dfAdult : df
          data frame of adult worms containing genotype information
     dfSel : df
          data fram of fitness benefit for each allele
     
     Returns
     ------
     selF : array
          array filling selF column, influencing fecundity
     selS : array      
          array filling selS column, influencing survival
     '''
     
     
def wormdf_fx(dfHost, villages, infhost, muWormBurden, sizeWormBurden, locus,
              initial_migration, initial_distance_m, theta, basepairs, mutation, 
              recombination, time2Ancestral, thetaRegional, time_join):
     '''builds worm dataframe
     Parameters
     ---------
     dfHost : df
          dataframe of host information
     villages : int
          number of villages, 2
     infhost : int,list
          number of infected hosts
     muWormBurden : int
          mu for neg binomial, 100
     sizeWormBurden : int
          size for neg binomial, 1
     locus : int
          how many loci to simulate
     initial_migration : float
          initial migration rat, 0.0001
     initial_distance_m : int
          distance between villages, [1000]
     theta : float, list
          theta value of diversity for loci and village, [[5,5],[10,10]]
     basepairs : int, list
          length of locus, [13000, 200000]
     mutation : float, list 
          mutation rate of each locus, [7.6E-8, 2.9E-9]
     recombination : float, list
          recombination rate of each locus, [0, 2.9E-9]
     time2ancestral : float
          generations to ancestral population, 1800
     thetaRegional : float
          diversity of regional population, 23
     time_join : float               
          time in past when 2 pops split, 240
          
     Returns
     -------
     dfAdult : df
          adult worms
     dfSel : df
          fitness contribution for each mutation

     '''
     dfAdult = pd.DataFrame({
                      'village' : [],
                      'hostidx' : [],
                      'age' : 1,
                      'sex' :  "F",
                      'R0net' : [],
                      'fec': [],
                      'loc0gt' : [],
                      'loc1gt' : [],
                      'selF' : [],
                      'selS' : []
                      })
     
     #create parasite burden per host
     popinit = []
     for mu, size, numworms in zip(muWormBurden, sizeWormBurden, infhost):
         #the number of worms per infected host
         wb_burden = np.random.negative_binomial(mu, mu/float(mu+size), numworms) 
         #[[10,12,1,15 ... infhostvill1],[19,1,5,4, ... infhostvill2]]
         popinit.append(np.array(wb_burden).tolist())
        
     #total worms per villages   
     wormpopsize = [sum(i) for i in popinit]
     
     #add village and hostidx to dfAdult
     dfvill = []
     dfhostidx = []
     for vill in range(villages):
          dfvill + (vill * wormpopsize[vill])
          for host in range(popinit[vill]):
               dfhostidx + ([dfHost[host]["hostidx"]] * popinit[vill][host])
     
     dfAdult["village"] = dfvill
     dfAdult["hostidx"] = dfhostidx
     
     #add age
     dfAdult["age"] = 1
     
     #add sex
     dfAdult["sex"] = random.choice("MF")
     
     #add R0net
     dfAdult["R0net"] = np.random.random(len(dfAdult.index))
     
     #add fec
     dfAdult["fec"] = 0
     
     #add genetic data
     posSel = []
     for loc in range(locus):
          gt_array, mutations = ms_outcall(wormpopsize, villages, initial_migration, 
                                    initial_distance_m, theta[loc], basepairs[loc], mutation[loc], 
                                    recombination[loc], time2Ancestral, thetaRegional, time_join)         
          posSel.append(mutations)
          dfAdult["locus_" + loc] = gt_array
     
     #create dfSel
     dfSel = sel_fx(posSel, 1)
          
     #add phenotype data
     selF, selS = fitness_fx(dfAdult, dfSel)
     dfAdult["selF"] = selF
     dfAdult["selS"] = selS        

     return dfAdult, dfSel
    

def wbsims_init(villages, villpopulation, prevalence, muTrans, sizeTrans, muWormBurden, 
                sizeWormBurden, locus, initial_migration, initial_distance_m, theta,
                basepairs, mutation, recombination, time2Ancestral, thetaRegional,
                time_join):
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
     infhost = villpopulation * prevalence
     #construct dfHost
     dfHost = host_fx(villages, infhost, muTrans, sizeTrans)
     #construct dfAdult, dfSel
     dfAdult, dfSel= wormdf_fx(dfHost, villages, infhost, muWormBurden, sizeWormBurden, 
                               locus, initial_migration, initial_distance_m, theta, 
                               basepairs, mutation, recombination, time2Ancestral, thetaRegional,
                               time_join)
     return dfAdult, dfHost, dfSel  
     

if __name__ == '__main__':
     #2 villages
     embed()
     wbsims_init(2, [100, 200], [0.1, 0.3], 100, 1, 5, 50, 2, 0.0001, [1000], 
                               [[5, 5], [10, 10]], [13000, 200000], 
                               [7.6E-8, 2.9E-9], [0, 2.9E-9], 1800, 23, 240)    
     
