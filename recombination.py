#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
'''
import numpy as np
import random
import pandas as pd

def recombination_fx(locus,
                     dfAdult,
                     recombination_rate,
                     basepairs):
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
              print("no males, no sex")
              break
         mf = 0
         while mf < dfAdult.loc[index, "fec"]:
              for loc in range(locus):
                   if recombination_rate[loc] == 0:
                        continue
                   else:
                        num_recomb = np.random.poisson(recombination_rate[loc] * basepairs[loc] * 2)
                        #print(num_recomb)
                        if num_recomb == 0:
                             row["locus_" + str(loc) + "_h1"] = row["locus_" + str(loc) + "_h" + random.choice("12")]
                             #male contribution
                             row["locus_" + str(loc) + "_h2"] = male["locus_" + str(loc) + "_h" + random.choice("12")]
                        else:
                             r = 0
                             #randomly choose male or female
                             sex_xing = random.choice("MF")
                             #while loop to account for multiple recombination events
                             h1m = male["locus_" + str(loc) + "_h1"]
                             h2m = male["locus_" + str(loc) + "_h2"]
                             h1f = row["locus_" + str(loc) + "_h1"]
                             h2f = row["locus_" + str(loc) + "_h2"]
                             print h1m
                             print h2m
                             print h1f
                             print h2f
                             if sex_xing is "M":
                                  while r < num_recomb:
                                       crossover_pos = random.randint(0, basepairs[loc])
                                       print crossover_pos
                                       try:
                                           hap1_co = [i for i,x in enumerate(h1m) if x > crossover_pos][-1]
                                           hap2_co = [i for i,x in enumerate(h2m) if x > crossover_pos][-1]
                                           h1m_new = h1m[0:hap1_co + 1] + h2m[hap2_co:]
                                           h2m_new = h2m[0:hap2_co + 1] + h1m[hap1_co:]
                                           h1m = h1m_new
                                           h2m = h2m_new
                                       except IndexError:
                                           continue
                                       r += 1
                             elif sex_xing is "F":
                                  while r < num_recomb:
                                       crossover_pos = random.randint(0, basepairs[loc])
                                       print crossover_pos
                                       try:
                                           hap1_co = [i for i,x in enumerate(h1f) if x > crossover_pos][-1]
                                           hap2_co = [i for i,x in enumerate(h2f) if x > crossover_pos][-1]
                                           h1f_new = h1f[0:hap1_co + 1] + h2f[hap2_co:]
                                           h2f_new = h2f[0:hap2_co + 1] + h1f[hap1_co:]
                                           h1f = h1f_new
                                           h2f = h2f_new
                                       except IndexError:
                                           continue
                                       r += 1
                             row["locus_" + str(loc) + "_h1"] = random.choice([h1f, h2f])
                             row["locus_" + str(loc) + "_h2"] = random.choice([h1m, h2m])
              dfAdult_mf = dfAdult_mf.append([row], ignore_index=True)
              mf += 1

    return dfAdult_mf
