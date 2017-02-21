#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""

import math
import subprocess
import numpy as np
import pandas as pd
import random
import pickle
import copy

from .agehost import agehost_fx
from .worm import Worms


def host_fx(village, infhost, muTrans, sizeTrans):
    '''Creates a transmission matrix for locations of infected hosts
       default_test : host_fx(2, [100, 300], 100, 1)

    Parameters
    ----------
    villages : int
        number of villages in the simulation
    infhost : list of pop. sizes
         number of infected hosts at intial time. prev*hostpopsize
    muTrans : float
        average distance in meters between hosts
    sizeTrans : float
        aggregation parameter

    Returns
    -------
    dfHost : dataframe
    '''
    deathdict = pickle.load( open( "../figs/data/acttable.p", "rb" ) )
    assert len(village) == len(infhost)
    coordinates = []
    host_idx = []

    for vill in range(len(village)):
         #list of host positions
         coordinates.extend( np.random.negative_binomial(sizeTrans, sizeTrans
                             / float((sizeTrans+muTrans)), (infhost[vill],
                                                           2)) + village[vill].dist)

         for host in range(infhost[vill]):
             host_idx.append("v" + str(vill) + "h" + str(host + 1))
    sex = [random.choice("01") for i in range(sum(infhost))]
    age_death = [agehost_fx(i, deathdict) for i in sex]

    dfHost = pd.DataFrame({
        'village': np.repeat(range(len(village)), infhost),
        'hostidx': host_idx,
        'sex': sex,
        'age': [i[0] for i in age_death],
        'agedeath': [i[1] for i in age_death],
        'coordinates': coordinates,
        'MDA': np.zeros(sum(infhost)),
        'MDA_cum': np.zeros(sum(infhost))
    })
    dfHost = dfHost.loc[:, ['village', 'hostidx', 'sex',
            'age', 'agedeath', 'coordinates', 'MDA', 'MDA_cum']]
    return(dfHost)


def coalsims_migmat_fx(numvillages, initial_migration, initial_distance_m, thetaN0,
                     basepairs, mutation_rate):
    '''Creates a string that represents a migration matrix between
    metapopulations. Migration here is a stepping stone not island model

    Format specified by ms (hudson 2000). Uses euclidian distances. The value
    of 2Nm is weighted by the distance from the next population as an
    exponential random variable.  The highest this can be is 2Nm/1

    Parameters
    ----------
    thetaN0 : float
         theta from the first village of the locus

    Returns
    -------
    Migration matrix : string
        migration matrix for scrm
    '''
    ne = thetaN0 / (4 * mutation_rate * basepairs)

    if numvillages > 4:
        raise ValueError("only handles 4 villages ATM")
    elif numvillages <= 4:
        assert len(initial_distance_m) == ((numvillages * (numvillages - 1)) / 2)
        mig = []  # initiate blank migration list
        for meters in initial_distance_m:
            mig.append((initial_migration) / (np.random.exponential(meters)))
        if numvillages == 2:
            m1 = 4 * ne * mig[0]  # 4Nm
            return("{}".format(m1))  # mig_matrix is symmetrical and island
        elif numvillages == 3:
            m1 = 4 * ne * mig[0]
            m2 = 4 * ne * mig[1]
            m3 = 4 * ne * mig[2]
            # mig_matrix is symmetrical
            return("{} {} {} {} {} {} {} {} {}".format(
                0, m1, m2, m1, 0, m3, m2, m3, 0))
        elif numvillages == 4:
            m1 = 4 * ne * mig[0]
            m2 = 4 * ne * mig[1]
            m3 = 4 * ne * mig[2]
            m4 = 4 * ne * mig[3]
            m5 = 4 * ne * mig[4]
            m6 = 4 * ne * mig[5]
            return("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(
                0, m1, m2, m3, m1, 0, m4, m5, m2, m4, 0, m6, m3, m5, m6, 0))


def parse_coalsims_fx(msout, ploidy, nind):
    """Parse output from ms or scrm.
    Parameters
    ---------
    msout : iterator,file
         stdout from scrm
    ploidy : int
         ploidy of the current locus

    Returns
    ------
    gt_array
         columns to be added to dfAdult
    pos
         list of positions for the mutaitons
    """
    print('NINDVS ' + str(nind))
    nind = int(nind)
    # Parsing positions
    for line in iter(msout.stdout.readline, ''):
        line = line.decode('utf-8')
        if line.startswith("positions"):
            # collisions can result here when theta is high
            pos = np.round(
                np.array(
                    line.strip().split()[1:],
                    dtype=np.float64))
            prev = 0
            # :TODO double check if this works correctly
            for idx, item in enumerate(pos, start=0):
                while prev >= item:
                    item += 1
                pos[idx] = item
                prev = pos[idx]
            break
        else:
            pass
    pos = pos.astype(np.int64)
    if ploidy == 1:
        gt_array = np.zeros((nind, pos.shape[0]) , dtype=np.uint8)
        #:TODO Seriously think about returning the whole genotype array
        cix = 0
        for line in iter(msout.stdout.readline, ''):
            line = line.decode('utf-8')
            line = list(line.strip())
            try:
                gt_array[cix, :] = np.array(line, dtype=np.uint8)
            except IndexError:
                break
            cix += 1
        return(gt_array, pos)
    elif ploidy == 2:
        cix = 0
        gt_array = np.zeros((nind, pos.shape[0]), dtype=np.uint8)
        gt_array2 = np.zeros((nind, pos.shape[0]) , dtype=np.uint8)
        for line in iter(msout.stdout.readline, ''):
            line = line.decode('utf-8')
            hap = np.array(list(line.strip()), dtype=np.uint8)
            gt_array[cix, :] = hap
            hap2 = next(iter(msout.stdout.readline, ''))
            gt_array2[cix, :] =  np.array(list(hap2.strip()),
                    dtype=np.uint8)
            cix += 1
            if cix == nind:
                break
        return(gt_array, gt_array2, pos)


def coalsims_fx(worm_popsize, numvillages, initial_migration, initial_distance_m,
        theta, basepairs, mutation_rate, recombination_rate, time2Ancestral, thetaRegional,
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
    gt2 : array
         genotype data for second haplotype
    mutations : array
         array of mutation positions
    '''
    # theta for first village
    thetaN0 = theta[0]
    #recombination rate for locus
    rho = thetaN0 * (recombination_rate / mutation_rate)
    rho is 0
    # time to the ancestral population
    tA = basepairs * (time2Ancestral / (thetaN0 / mutation_rate))
    # time to join villages
    tjoin = basepairs * (time_join / (thetaN0 / mutation_rate))
    #check ploidy for total haplotypes
    if rho == 0:
        ploidy = 1
    else:
        worm_popsize = [x*2 for x in worm_popsize]
        ploidy = 2
    # parameters for ms or scrm
    ms_params = {
        'nhaps': sum(worm_popsize),
        'theta': thetaN0,
        'rho': rho,
        'basepair': basepairs - 1,
        'exp_growth': (-1 / tA) * math.log(thetaRegional),
        'time_growth': tA,
        'sig_digits': 12,
        'tjoin' : tjoin,
        'seed' : 12345
    }

    # order matters for the command call
    if numvillages == 1:
        scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                     "-eG {join} {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                     " {sig_digits} -seed {seed} ")
        mscmd = scrm_base.format(**ms_params)
    else:  # ms setup for >1 villages
        num_subpops = len(worm_popsize)  # -I num_pops
        sub_pop = " ".join(map(str, worm_popsize))  # -I X i j ...
        mm = coalsims_migmat_fx(
           numvillages,
           initial_migration,
           initial_distance_m,
           thetaN0,
           basepairs,
           mutation_rate)
        if numvillages == 2:
            scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                         "{sub_pop} {present_pop} {join} -G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                         " {sig_digits} ")
            ms_params[
                'present_pop'] = '-n 1 1 -n 2 {}'.format(float(theta[1]) / thetaN0)
            ms_params[
                'sub_pop'] = '-I {} {} {}'.format(num_subpops, sub_pop, mm)
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

            for village_ix in range(numvillages):
                present_pop = float(theta[village_ix]) / thetaN0
                subpopstr += '-n {} {} '.format(village_ix + 1, present_pop)
                if village_ix != numvillages - 1:
                    joinstr += '-ej {0} {1} {2} '.format(
                        tjoin,
                        village_ix + 1,
                        village_ix + 2)
            ms_params['present_pop'] = subpopstr
            ms_params['join'] = joinstr
            mscmd = scrm_base.format(**ms_params)
    print(mscmd)
    msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
    #msout = subprocess.check_output(mscmd.split(" ")).decode('utf-8')
    #:TODO refactor this, less variation on returns
    if ploidy == 1:
        gt, positions = parse_coalsims_fx(msout, ploidy, sum(worm_popsize))
        return(gt, positions)
    elif ploidy == 2:
        gt, gt2, positions = parse_coalsims_fx(msout, ploidy,
                sum(worm_popsize)/2)
        return(gt, gt2, positions)


def sel_fx(dfAdult, basepairs, perc_locus, cds_length, intgen_length):
    '''Initializes the distribution of fitness effects for each mutation

    Parameters
    ---------
    positions : list
      list from coalsims_fx, the the line "positions" in the scrm/ms output

    Returns
    ------
    dfSel : df
        df containing mutations and fitness effect
    cds_coordinates : nested array
        coordinates (start, end) of coding sequences

    '''
    #size parameter
    size = 3
    #average distance between in bp
    mu = intgen_length
    #last coordinate is approx num_cds * mu; so if num_cds is too large or mu is too long
    #genes will run over the locus length
    #dfAdult.pos
    for locus in dfAdult.h2.keys():
  