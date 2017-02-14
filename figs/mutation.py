#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:59:27 2017

@author: scott
"""
import numpy as np
import random

def mutation_fx(locus,
                dfAdult_mf,
                mutation_rate,
                recombination_rate,
                basepairs):
    '''Calculates number of mutations, changes haplotypes

    Parameters
    ---------
    locus: int
         number of loci
    dfAdult_mf : figs.Worms object 
          New larval parasites
    mutation_rate : float, list
          mutation rates for each locus
    recombination_rate : float, list
        recombination rates for each locus
    basepairs : int, list
          length of each locus in basepairs
    cds_coordinates : list
        list of coding positions

    Returns
    ------
    dfAdult_mf : fig.Worms object
         larval parasites with mutated positions
    mutations : int, list
         list of positions of new mutations
    '''
    #dfAdult_mf.reset_index(drop=True, inplace=True)
    new_positions = []
    nworms = dfAdult_mf.meta.shape[0]
    for loc in range(locus):
        if recombination_rate[loc] == 0:
            mut_coef = 1
        else: 
            mut_coef = 2
        num_muts = np.random.binomial(mut_coef * nworms, 
                basepairs[loc] * mutation_rate[loc])
        #print num_muts
        positions = dfAdult_mf.pos["locus_" + loc]
        for mut in range(num_muts):
            randmf = np.random.randint(0, nworms)
            narray = np.zeros(nworms, np.uint8)
            narray[randmf] = 1
            newsite = np.random.randint(1, basepairs[loc])
            iix = np.argmax(positions > newsite)
            if recombination_rate[loc] == 0:
                dfAdult_mf.h1["locus_" + str(loc)] =\
                        np.insert(dfAdult_mf.h1["locus_" + str(loc)], 
                                iix, narray, axis=1)
            else:
                oarray = np.zeros(nworms, np.uint8)
                whap = random.choice([1,2])
                if whap == 1: whap2 = 2
                else: whap2 = 1
                hap =\
                        getattr(dfAdult_mf,"h"+\
                        str(whap))["locus_"+str(loc)]  
                hap = np.insert(hap, iix, narray, axis=1)
                ohap =\
                        getattr(dfAdult_mf, "h"+\
                        str(whap2))["locus_"+str(loc)]
                ohap = np.insert(ohap, iix, oarray, axis=1)
            new_positions.append(newsite)
    return(dfAdult_mf, new_positions)







