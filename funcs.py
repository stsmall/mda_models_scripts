import time
import math
import subprocess
from collections import defaultdict
import re
import numpy as np
import random


def migration_matrix(villages, initial_migration, initial_distance_m, theta, 
        basepairs, mutation):
    '''Creates a string that represents a migration matrix between
    metapopulations.

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

    Returns
    -------
    Migration matrix: string
        migration matrix for scrm
    '''

    print("starting migration_matrix")
    ne = theta / (4*mutation*basepairs)
    #cant figure out how to pythonically increase the villages without just using loops
    if villages > 4:         
        raise ValueError("only handles 4 villages ATM")
    elif villages < 4:
        if len(initial_distance_m) != ((villages)*(villages-1)/2): 
            raise ValueError(("there are not adequate pairwise comparisons in" 
            "distance_m to match villages"))
        mig = [] # initiate blank migration list
        for meters in initial_distance_m:
            mig.append((initial_migration)/(np.random.exponential(meters)))
        if villages == 2:
            m1 = 4*ne*mig[0] #4Nm
            return "{}".format(m1) #mig_matrix is symmetrical and island
        elif villages == 3:
            m1 = 4*ne*mig[0]
            m2 = 4*ne*mig[1]
            m3 = 4*ne*mig[2]
            return "{} {} {} {} {} {} {} {} {}".format(0, m1, m2, m1, 
                    0, m3, m2, m3, 0) #mig_matrix is symmetrical
        elif villages == 4:
            m1 = 4*ne*mig[0]
            m2 = 4*ne*mig[1]
            m3 = 4*ne*mig[2]
            m4 = 4*ne*mig[3]
            m5 = 4*ne*mig[4]
            m6 = 4*ne*mig[5]
            return "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(0, m1, 
                    m2, m3, m1, 0, m4, m5, m2, m4, 0, m6, m3, m5, m6, 0)


def parse_ms_output():
    """Parse output from ms or scrm.
    """
    pass


def ms_outcall(worm_popsize, villages, initial_migration, initial_distance_m, 
        theta, basepairs, mutation, recombination, time2Ancestral, 
        thetaAncestral,  thetaRegional, time_join12, time_join23, time_join34):
    '''External call to ms (Hudson 2000) or scrm. 
    
    Calls the function migration_matrix() if using more than 1 village. 
    This function will then read in the stdout from ms/scrm. The 
    location of segsites are assigned a random number then scaled
    by the length of the desired sequence in basepairs.

    Parameters 
    ----------
    worm_popsize: (list,int) 
        how many worms to simulate
    villages: int 
        number of villages/metapopulations
    theta: (list,float) 
        list of theta values for each locus. With more than 1 village 
        it is list of lists [[locus1_meta1,locus1_meta2],[locus2_meta1,locus2_meta2]]
    basepairs:(list,int) 
        list of basepairs (lengths) for each locus
    mutation: (list,float) 
        mutation rate for each locus as probability per base per generation
    recombination:(list,float) 
        recombination rate for each locus as probability per base per generation
    thetaAncestral:float 
        theta for ancestral pops to the ratio of N0; e.g.,
    thetaRegional: float 
        theta for regional pops to the ratio of N0 e.g.,23 times larger
    time2Ancestral: int 
        time in generations to the ancestral population
    time_join12: int 
        time in generations for joining/splitting pop 1 into 2
    time_join23: int 
        time in generations for joining/splitting pop 2 into 3
    time_join34: int 
        time in generations for joining/splitting pop 3 into 4

    Returns
    -------


    '''
    print("starting scrm...")
    t0=time.clock()
    num_loci = len(theta) 
    #theta for first population, all other thetas are scaled from this value
    thetaN0 = [t[0] for t in theta]     
    #population recombination rate
    rho = [t*(r/u) for t, r, u in zip(thetaN0, recombination, mutation)] 
    # :TODO simplify this
    #time to the ancestral population
    tA = [b*(time2Ancestral/(t/m)) for b, t, m in zip(basepairs, thetaN0, mutation)] 
    # time to join 1 & 2
    t12 = [b*(time_join12/(t/m)) for b, t, m in zip(basepairs, thetaN0, mutation)]     
    # time to join 2 & 3
    t23 = [b*(time_join23/(t/m)) for b, t, m in zip(basepairs, thetaN0, mutation)]     
    # time to join 3 & 4
    t34 = [b*(time_join34/(t/m)) for b, t, m in zip(basepairs, thetaN0, mutation)]
    #list intialization for recording haplotypes
    hap_pop = defaultdict(list) 
    time_join = [t12, t23, t34]

    for loc in range(0, num_loci):
        if rho[loc] is 0:
            ploidy = 1
        else:
            ploidy = 2
        worm_popsize[:] = [x * ploidy for x in worm_popsize]
        total_inds = sum(worm_popsize) # how many
        #positions = '%.2E' %(Decimal(str(basepairs[i])))
        ms_params = {
                'nhaps': total_inds,
                'theta': thetaN0[loc], 
                'rho': rho[loc],
                'basepair': basepairs[loc] - 1,
                'exp_growth': (-1/tA[loc])*math.log(thetaRegional),
                'time_growth': tA[loc],
                'sig_digits': 12,
                }
        scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} " 
                     "-G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                     " {sig_digits}")
        
        if villages == 1: 
            mscmd = scrm_base.format(**ms_params)
        else: #ms setup for >1 villages
            num_subpops = len(worm_popsize) #-I num_pops
            sub_pop = " ".join(map(str, worm_popsize))#-I X i j ...
            mscmd = scrm_base
            for village_ix in range(villages):
                mm = migration_matrix(villages, initial_migration,
                    initial_distance_m, thetaN0[loc], basepairs[loc],
                    mutation[loc])
                if village_ix == 0: 
                    present_pop = 1
                else: 
                    present_pop = float(thetaN0[loc])/theta[loc][village_ix]
                    mscmd += '-ej {0} {1} {2}'.format(time_join[village_ix-1][loc], 
                            village_ix + 1, village_ix + 2)
                mscmd += '-n {} {} '.format(village_ix, present_pop) 
            if villages == 2:
                # Island model
                mscmd += '-I {} {} {} '.format(num_subpops, sub_pop, mm)
            else:
                mscmd += 'I {} {} '.format(num_subpops, sub_pop)
                mscmd += '-ma {} '.format(mm)
        print(mscmd)
        proc = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
        #parses ms output from stdout
        for line in iter(proc.stdout.readline, ''):
            if line.startswith("positions"):
                positions = [round(i) for i in map(float, line.strip().split()[1:])]
                for line in iter(proc.stdout.readline, ''):
                    hap = line.strip()
                    hap_pop["locus"+str(loc)].append([positions[j] for j in\
                            [m.start() for m in re.finditer("1", hap)]])  
    print(time.clock()-t0)
    return hap_pop
    #{'locus0':[14.0, 26.0],[5.0, 7.0]]}


def recombination_fx(mf, num_recomb, basepairs):
    """this is run every time the prob of a recombination is true

    Parameters
    ----------
    mf: list
        list of mf from adults
    num_recomb: list of ints
        number of recombination events observed

    Returns
    -------
    """
    for i, bp in enumerate(basepairs): #number of muts
        recomb = 0
        while recomb < num_recomb[i]: #keep going until all recombinations are  assigned
            rec_mf = random.randrange(0, len(mf)) #choose random index in mf lis
            new_recomb = random.randint(0, 3)
            if new_recomb < 2: #first parent
                hap1 = mf[rec_mf][i+1][0]
                hap2 = mf[rec_mf][i+1][1]
                crossover_pos = random.randint(0, bp)
                hap1.sort()
                hap2.sort()
                hap1_co = next(l[0] for l in enumerate(hap1) if l[1] > crossover_pos)
                hap2_co = next(l[0] for l in enumerate(hap2) if l[1] > crossover_pos)
                hap1_new = hap1[0:hap1_co] + hap2[hap2_co:]
                hap2_new = hap2[0:hap2_co] + hap1[hap1_co:]
                mf[rec_mf][i+1][0] = hap1_new
                mf[rec_mf][i+1][1] = hap2_new
            else:#second parent
                hap3 = mf[rec_mf][i+1][2]
                hap4 = mf[rec_mf][i+1][3]
                crossover_pos = random.randint(0, bp)
                hap3.sort()
                hap4.sort()
                hap3_co = next(l[0] for l in enumerate(hap3) if l[1] > crossover_pos)
                hap4_co = next(l[0] for l in enumerate(hap4) if l[1] > crossover_pos)
                hap3_new = hap1[0:hap3_co] + hap2[hap4_co:]
                hap4_new = hap2[0:hap4_co] + hap1[hap3_co:]
                mf[rec_mf][i+1][2] = hap3_new
                mf[rec_mf][i+1][3] = hap4_new
            recomb += 1
    return mf
