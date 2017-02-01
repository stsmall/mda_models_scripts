#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:59:27 2017

@author: scott
"""
import numpy as np
import random
import pandas as pd

def mutation_fx(locus, dfAdult_mf, mutation_rate, recombination_rate, basepairs, cds_coordinates):
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
              muts = 0     
              if num_muts != 0:
                   while muts < num_muts:
                        randmf = np.random.randint(0,len(dfAdult_mf))
                        newsite = np.random.randint(1,basepairs[loc])
                        #positions.append(newsite)
                        newhap = np.append(dfAdult_mf.iloc[randmf]["locus_" + str(loc)], newsite)
                        dfAdult_mf.set_value(randmf, "locus_" + str(loc), newhap) 
                        muts += 1
          
         else:
              num_muts = np.random.binomial(2 * len(dfAdult_mf), basepairs[loc] * mutation_rate[loc])
              #print num_muts
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
                   cds_positions = []
                   for start, end in cds_coordinates[loc - 1]:
                        cds_positions.extend(positions[np.where(np.logical_and(positions >= start, 
                                            positions <= end))])                      
                   muts_counter.append(cds_positions) 
    dfMuts = pd.DataFrame({
                           "locus" : np.repeat(range(1, locus), muts_counter),
                           "position" : sum(muts_counter, []),
                           "selF" :  np.zeros(len(positions)),
                           "selS" :  np.zeros(len(positions)),
                           "freqInit" :  np.zeros(len(positions))
                           })
    return dfAdult_mf, dfMuts                    
                        
                        
                        
                        
                        
                        
                        
                        