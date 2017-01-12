import math
import subprocess
from collections import defaultdict
import re
import numpy as np


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

    Returns
    -------
    Migration matrix: string
        migration matrix for scrm
    '''
    ne = theta / (4 * mutation * basepairs)
    # cant figure out how to pythonically increase the villages without just
    # using loops
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


def parse_ms_output():
    """Parse output from ms or scrm.
    """
    pass


def ms_outcall(
        worm_popsize,
        villages,
        initial_migration,
        initial_distance_m,
        theta,
        basepairs,
        mutation,
        recombination,
        time2Ancestral,
        thetaRegional,
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
    '''
    num_loci = len(theta)
    # theta for first population, all other thetas are scaled from this value
    thetaN0 = [t[0] for t in theta]
    # population recombination rate
    rho = [t * (r / u) for t, r, u in zip(thetaN0, recombination, mutation)]
    # :TODO simplify this
    # time to the ancestral population
    tA = [b * (time2Ancestral / (t / m))
          for b, t, m in zip(basepairs, thetaN0, mutation)]
    # time to join villages
    tjoin = [b * (time_join / (t / m))
           for b, t, m in zip(basepairs, thetaN0, mutation)]
    # list intialization for recording haplotypes
    hap_pop = defaultdict(list)

    for loc in range(0, num_loci):
        if rho[loc] is 0:
            ploidy = 1
        else:
            ploidy = 2
        worm_popsize[:] = [x * ploidy for x in worm_popsize]
        total_inds = sum(worm_popsize)  # how many
        #positions = '%.2E' %(Decimal(str(basepairs[i])))
        ms_params = {
            'nhaps': total_inds,
            'theta': thetaN0[loc],
            'rho': rho[loc],
            'basepair': basepairs[loc] - 1,
            'exp_growth': (-1 / tA[loc]) * math.log(thetaRegional),
            'time_growth': tA[loc],
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
               thetaN0[loc],
               basepairs[loc],
               mutation[loc])
            
            if villages == 2:
                 scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                     "{sub_pop} {present_pop} {join} -G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                     " {sig_digits} ")
                 ms_params['present_pop'] = '-n 1 1 -n 2 {}'.format(float(theta[loc][1])/thetaN0[loc])
                 ms_params['sub_pop'] = '-I {} {} {}'.format(num_subpops, sub_pop, mm)
                 ms_params['join'] = '-ej {} 1 2'.format(tjoin[loc])
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
                    present_pop = float(theta[loc][village_ix])/thetaN0[loc]
                    subpopstr += '-n {} {} '.format(village_ix + 1, present_pop)
###dont want last set of -ej time 3 4 # for villages of 3 since there is no 4 to join
                    if village_ix != villages -1:
                         joinstr += '-ej {0} {1} {2} '.format(
                             tjoin[loc],
                             village_ix + 1,
                             village_ix + 2)
###this was a sloppy way of avoiding it, any ideas?
                 ms_params['present_pop'] = subpopstr
                 ms_params['join'] = joinstr
                 mscmd = scrm_base.format(**ms_params)     
        
        print(mscmd)
        proc = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)

        # parses ms output from stdout
        for line in iter(proc.stdout.readline, ''):
            if line.startswith("positions"):
                positions = [
                    round(i) for i in map(
                        float,
                        line.strip().split()[
                            1:])]
                for line in iter(proc.stdout.readline, ''):
                    hap = line.strip()
                    ti = [positions[j] for j in [m.start() 
                        for m in re.finditer("1", hap)]]
                    hap_pop["locus" + str(loc)].append(ti)
    return hap_pop
    #{'locus0':[14.0, 26.0],[5.0, 7.0]]}

if __name__ == '__main__':
#2villages
     ms_outcall([100,100], 2, 0.0001, [1000], [[5, 5], [10, 10]], [13000, 200000], 
                [7.6E-8, 2.9E-9],  [0, 2.9E-9], 1800, 23, 240)
#3villages
     #ms_outcall([100,100,100], 3, 0.0001, [1000,2000,500], [[5, 2, 2], 
      #           [10, 12, 12]], [13000, 200000], [7.6E-8, 2.9E-9],  
       #          [0, 2.9E-9], 1800, 23, 240)