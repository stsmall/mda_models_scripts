#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:59:27 2017

@author: scott
"""
import numpy as np
import random
import pandas as pd

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
    position = []
    muts_counter = []      
    for loc in range(locus):
         if recombination_rate[loc] == 0:
              num_muts = np.random.binomial(len(dfAdult_mf), basepairs[loc] * mutation_rate[loc])
              muts_counter.append(num_muts)
              muts = 0     
              if num_muts != 0:
                   while muts < num_muts:
                        randmf = np.random.randint(0,len(dfAdult_mf))
                        newsite = np.random.randint(1,basepairs[loc])
                        position.append(newsite)
                        newhap = np.append(dfAdult_mf.iloc[randmf].locus, newsite)
                        dfAdult_mf.set_value(randmf, "locus_" + str(loc), newhap) 
                        muts += 1
          
         else:
              num_muts = np.random.binomial(2 * len(dfAdult_mf), basepairs[loc] * mutation_rate[loc])
              muts_counter.append(num_muts)
              muts = 0     
              if num_muts != 0:
                   while muts < num_muts:
                        randmf = np.random.randint(0,len(dfAdult_mf))
                        newsite = np.random.randint(1,basepairs[loc])
                        position.append(newsite)
                        newhap = np.append(dfAdult_mf.iloc[randmf].locus, newsite)
                        dfAdult_mf.set_value(randmf, "locus_" + str(loc) + "_h" + random.choice("12"), newhap) 
                        muts += 1
                        
    dfMuts = pd.DataFrame({
                           "locus" : np.repeat(range(locus), muts_counter),
                           "position" : position,
                           "selF" : 0,
                           "selS" : 0,
                           "freqInit" : 0
                           })
    return dfAdult_mf, dfMuts                    
                        
                        
                        
                        
                        
                        
                        
                        