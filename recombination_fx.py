import numpy as np
import random
import pandas as pd

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
         male = dfAdult[dfAdult.sex == "M"].sample(1) #choose a random male
         mf = 0
         while mf < dfAdult.loc[index, "fec"]:
              for loc in range(locus):
                   if recombination_rate[loc] == 0:
                        pass
                   else:                                         
                        num_recomb = np.random.poisson(recombination_rate[loc] * basepairs[loc] * 2)
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
                             if sex_xing is "M":
                                  while r < num_recomb:
                                       crossover_pos = random.randint(0, basepairs[loc])
                                       hap1_co = next(l[0] for l in enumerate(
                                                 h1m) if l[1] > crossover_pos)
                                       hap2_co = next(l[0] for l in enumerate(
                                                 h2m) if l[1] > crossover_pos)
                                       h1m_new = h1m[0:hap1_co] + h2m[hap2_co:] 
                                       h2m_new = h2m[0:hap2_co] + h1m[hap1_co:]
                                       h1m = h1m_new
                                       h2m = h2m_new
                                       r += 1
                             elif sex_xing is "F":
                                  while r < num_recomb:
                                       crossover_pos = random.randint(0, basepairs[loc])
                                       hap1_co = next(l[0] for l in enumerate(
                                                 h1f) if l[1] > crossover_pos)
                                       hap2_co = next(l[0] for l in enumerate(
                                                 h2f) if l[1] > crossover_pos)
                                       h1f_new = h1f[0:hap1_co] + h2f[hap2_co:] 
                                       h2f_new = h2f[0:hap2_co] + h1f[hap1_co:]
                                       h1f = h1f_new
                                       h2f = h2f_new
                                       r += 1
                             row["locus_" + str(loc) + "_h1"] = random.choice([h1f, h2f])     
                             row["locus_" + str(loc) + "_h2"] = random.choice([h1m, h2m])      
              dfAdult_mf = dfAdult_mf.append([row], ignore_index=True)
              mf += 1
     
    return dfAdult_mf

#                crossover_pos = random.randint(0, bp)
#                hap1_co = next(l[0] for l in enumerate(
#                    hap1) if l[1] > crossover_pos)
#                hap2_co = next(l[0] for l in enumerate(
#                    hap2) if l[1] > crossover_pos)
#                hap1_new = hap1[0:hap1_co] + hap2[hap2_co:]
#                hap2_new = hap2[0:hap2_co] + hap1[hap1_co:]