import random
import math

import numpy as np


def hostdeath(meta_popdict, mpop):
    '''when a host dies/cures all the MF/Wb disappear from the dict'''
    dpop = random.choice(meta_popdict[mpop].keys()) #draw random host infection
    meta_popdict[mpop][dpop] = {} #blanks the dictionary entry
    return meta_popdict


def maturation(mpop, meta_popdict, time_month, density_dependence, 
        mortalityHost, mortalityAdult, mortalityJuv, mortalityMF, 
        fecundity, basepairs, mutation, recombination):
    """Life cycle stage, including births/deaths of worms

    Parameters
    ----------
    mpop: int
        metapopulation counter
    meta_popdict: dict 
        perform functions on dictionary containing vill,host,worm and genotypes
    time_month: int 
        what month is it 1 ... 12
    infhost: int 
        hostpopsize * prev = number of infected hosts
    density_dependence: boolean 
        if True then use density dependencec functions to calculate mortality/fecundity
    mortalityHost: float 
        prob of host dying per year, default 1/70
    mortalityAdult: float 
        prob of adult surviving per year
    mortalityJuv: float 
        prob of Juv surviving per month
    mortalityMF: float 
        prob of MF surviving per month
    fecundity: int 
        average number of MF produced per adult per month

    Returns
    -------
    """
    mpop = "meta_" + str(mpop+1)
    print(mpop)
    mf_sum = []
    #since this month to month
    if time_month%12 is 0: #this denotes 1 year has passed so adults mature to next age class
        if random.uniform(0, 1) < mortalityHost:
            meta_popdict = hostdeath(meta_popdict, mpop)
        for npop in meta_popdict[mpop].keys(): #infected hosts as pops within villages
            #count individuals in each class for density dependent calculations
            sum_adult = sum([len(meta_popdict[mpop][npop][item]) for item in ['A_1', 'A_2', 'A_3', 'A_4', 'A_5', 'A_6', 'A_7', 'A_8']])
            sum_juv = sum([len(meta_popdict[mpop][npop][item]) for item in ['J_1', 'J_2', 'J_3', 'J_4', 'J_5', 'J_6', 'J_7', 'J_8', 'J_9', 'J_10', 'J_11', 'J_12']])
            sum_mf = sum([len(meta_popdict[mpop][npop][item]) for item in ['MF_1', 'MF_2', 'MF_3', 'MF_4', 'MF_5', 'MF_6', 'MF_7', 'MF_8', 'MF_9', 'MF_10', 'MF_11', 'MF_12']])
            mf_sum.append(sum_mf)

            if (sum_adult == 0) and (sum_juv == 0) and (sum_mf == 0):
                meta_popdict[mpop][npop] = {}
                #break

            if density_dependence: # density dependent mortality
                Km = 100  #carrying capacity
                bm = 1.0/Km #Km is carrying capacity
                Sm = 0.99 #maximum survival
                am = (Sm*bm)/math.exp(-1)
                mort_A = am * sum_adult * math.exp(-bm * sum_adult) #Ricker fx
                mort_J = am * sum_juv * math.exp(-bm * sum_juv) #Ricker fx
                mort_M = am * sum_mf * math.exp(-bm * sum_mf) #Ricker fx
                fecund = fecundity #fecundity per month
            else: # fixed mortality, actually prob of surviving
                mort_A = mortalityAdult #year
                mort_J = mortalityJuv #month
                mort_M = mortalityMF #month
                fecund = fecundity #fecundity per month

           # adult age classes
            meta_popdict[mpop][npop]["A_8"] = random.sample(meta_popdict[mpop][npop]["A_7"], int(round(len(meta_popdict[mpop][npop]["A_7"])*mort_A)))
            meta_popdict[mpop][npop]["A_7"] = random.sample(meta_popdict[mpop][npop]["A_6"], int(round(len(meta_popdict[mpop][npop]["A_6"])*mort_A)))
            meta_popdict[mpop][npop]["A_6"] = random.sample(meta_popdict[mpop][npop]["A_5"], int(round(len(meta_popdict[mpop][npop]["A_5"])*mort_A)))
            meta_popdict[mpop][npop]["A_5"] = random.sample(meta_popdict[mpop][npop]["A_4"], int(round(len(meta_popdict[mpop][npop]["A_4"])*mort_A)))
            meta_popdict[mpop][npop]["A_4"] = random.sample(meta_popdict[mpop][npop]["A_3"], int(round(len(meta_popdict[mpop][npop]["A_3"])*mort_A)))
            meta_popdict[mpop][npop]["A_3"] = random.sample(meta_popdict[mpop][npop]["A_2"], int(round(len(meta_popdict[mpop][npop]["A_2"])*mort_A)))
            meta_popdict[mpop][npop]["A_2"] = random.sample(meta_popdict[mpop][npop]["A_1"], int(round(len(meta_popdict[mpop][npop]["A_1"])*mort_A)))
            # first year adults
            meta_popdict[mpop][npop]["A_1"] = random.sample(meta_popdict[mpop][npop]["J_12"], int(round(len(meta_popdict[mpop][npop]["J_12"])*mort_J)))
            # juvenille
            meta_popdict[mpop][npop]["J_12"] = random.sample(meta_popdict[mpop][npop]["J_11"], int(round(len(meta_popdict[mpop][npop]["J_11"])*mort_J)))
            meta_popdict[mpop][npop]["J_11"] = random.sample(meta_popdict[mpop][npop]["J_10"], int(round(len(meta_popdict[mpop][npop]["J_10"])*mort_J)))
            meta_popdict[mpop][npop]["J_10"] = random.sample(meta_popdict[mpop][npop]["J_9"], int(round(len(meta_popdict[mpop][npop]["J_9"])*mort_J)))
            meta_popdict[mpop][npop]["J_9"] = random.sample(meta_popdict[mpop][npop]["J_8"], int(round(len(meta_popdict[mpop][npop]["J_8"])*mort_J)))
            meta_popdict[mpop][npop]["J_8"] = random.sample(meta_popdict[mpop][npop]["J_7"], int(round(len(meta_popdict[mpop][npop]["J_7"])*mort_J)))
            meta_popdict[mpop][npop]["J_7"] = random.sample(meta_popdict[mpop][npop]["J_6"], int(round(len(meta_popdict[mpop][npop]["J_6"])*mort_J)))
            meta_popdict[mpop][npop]["J_6"] = random.sample(meta_popdict[mpop][npop]["J_5"], int(round(len(meta_popdict[mpop][npop]["J_5"])*mort_J)))
            meta_popdict[mpop][npop]["J_5"] = random.sample(meta_popdict[mpop][npop]["J_4"], int(round(len(meta_popdict[mpop][npop]["J_4"])*mort_J)))
            meta_popdict[mpop][npop]["J_4"] = random.sample(meta_popdict[mpop][npop]["J_3"], int(round(len(meta_popdict[mpop][npop]["J_3"])*mort_J)))
            meta_popdict[mpop][npop]["J_3"] = random.sample(meta_popdict[mpop][npop]["J_2"], int(round(len(meta_popdict[mpop][npop]["J_2"])*mort_J)))
            meta_popdict[mpop][npop]["J_2"] = random.sample(meta_popdict[mpop][npop]["J_1"], int(round(len(meta_popdict[mpop][npop]["J_1"])*mort_J)))
            # microfilaria
            meta_popdict[mpop][npop]["MF_12"] = random.sample(meta_popdict[mpop][npop]["MF_11"], int(round(len(meta_popdict[mpop][npop]["MF_11"])*mort_M)))
            meta_popdict[mpop][npop]["MF_11"] = random.sample(meta_popdict[mpop][npop]["MF_10"], int(round(len(meta_popdict[mpop][npop]["MF_10"])*mort_M)))
            meta_popdict[mpop][npop]["MF_10"] = random.sample(meta_popdict[mpop][npop]["MF_9"], int(round(len(meta_popdict[mpop][npop]["MF_9"])*mort_M)))
            meta_popdict[mpop][npop]["MF_9"] = random.sample(meta_popdict[mpop][npop]["MF_8"], int(round(len(meta_popdict[mpop][npop]["MF_8"])*mort_M)))
            meta_popdict[mpop][npop]["MF_8"] = random.sample(meta_popdict[mpop][npop]["MF_7"], int(round(len(meta_popdict[mpop][npop]["MF_7"])*mort_M)))
            meta_popdict[mpop][npop]["MF_7"] = random.sample(meta_popdict[mpop][npop]["MF_6"], int(round(len(meta_popdict[mpop][npop]["MF_6"])*mort_M)))
            meta_popdict[mpop][npop]["MF_6"] = random.sample(meta_popdict[mpop][npop]["MF_5"], int(round(len(meta_popdict[mpop][npop]["MF_5"])*mort_M)))
            meta_popdict[mpop][npop]["MF_5"] = random.sample(meta_popdict[mpop][npop]["MF_4"], int(round(len(meta_popdict[mpop][npop]["MF_4"])*mort_M)))
            meta_popdict[mpop][npop]["MF_4"] = random.sample(meta_popdict[mpop][npop]["MF_3"], int(round(len(meta_popdict[mpop][npop]["MF_3"])*mort_M)))
            meta_popdict[mpop][npop]["MF_3"] = random.sample(meta_popdict[mpop][npop]["MF_2"], int(round(len(meta_popdict[mpop][npop]["MF_2"])*mort_M)))
            meta_popdict[mpop][npop]["MF_2"] = random.sample(meta_popdict[mpop][npop]["MF_1"], int(round(len(meta_popdict[mpop][npop]["MF_1"])*mort_M)))

            #add generation to A_1 since they are new from Juv12 and 1 year has passed so all are moving to A2
            for subl in meta_popdict[mpop][npop]["A_1"]:
                subl[0] += 1

            # birth of MF_1
            mf1 = []
            for i in range(1, 8):
                for ad in meta_popdict[mpop][npop]["A_{}".format(i)]:
                    births = np.random.poisson(fecund)
                    n = 0
                    age_class = random.randint(1,8)
                    while n < births:
                        newmf1 = copy.copy(ad)
                        wb_parent2 = meta_popdict[mpop][npop]["A_{}".format(age_class)][random.randint(0,len(meta_popdict[mpop][npop]["A_{}".format(age_class)]))] #random adult
                        while len(wb_parent2) is 0:
                             age_class = random.randint(1,8)  
                             meta_popdict[mpop][npop]["A_{}".format(age_class)][random.randint(0,len(meta_popdict[mpop][npop]["A_{}".format(age_class)]))] #random adult
                        for p in range(2, len(newmf1)):
                            newmf1[p] = newmf1[p] + wb_parent2[p]
                        mf1.append(newmf1)
                        n += 1
            mf = mf1
           #mf = sum(mf1, [])

           # recombination in new mf
            num_recomb = []
            for i, bp in enumerate(basepairs):
                num_recomb.append(np.random.binomial(2*len(mf), (bp-1) * recombination[i]))
            if sum(num_recomb) != 0:
                mf = recombination_fx(mf, num_recomb, basepairs)

           #makes diploid
            for m in mf:
                for item in m:
                    if len(item) == 4:
                        item = [item[random.randint(0, 1)], item[random.randint(2, 3)]]

            #mutation in new mf
            num_muts = []
            for i in range(0, len(basepairs)):
                if recombination[i] == 0:
                    num_muts.append(np.random.binomial(len(mf), basepairs[i] * mutation[i]))
                else:
                    num_muts.append(np.random.binomial(2*len(mf), basepairs[i] * mutation[i]))
            if sum(num_muts) != 0:
                mf = mutation_fx(mf, num_muts, basepairs)
            
            meta_popdict[mpop][npop]["MF_1"] = mf

    else: #a year has not passed on months, juveniles and MF move to next age class
        for npop in meta_popdict[mpop].keys(): #inf
            print(npop)
            #count individuals in each class for density dependent calculations
            sum_adult = sum([len(meta_popdict[mpop][npop][item]) for item in ['A_1', 'A_2', 'A_3', 'A_4', 'A_5', 'A_6', 'A_7', 'A_8']])
            sum_juv = sum([len(meta_popdict[mpop][npop][item]) for item in ['J_1', 'J_2', 'J_3', 'J_4', 'J_5', 'J_6', 'J_7', 'J_8', 'J_9', 'J_10', 'J_11', 'J_12']])
            sum_mf = sum([len(meta_popdict[mpop][npop][item]) for item in ['MF_1', 'MF_2', 'MF_3', 'MF_4', 'MF_5', 'MF_6', 'MF_7', 'MF_8', 'MF_9', 'MF_10', 'MF_11', 'MF_12']])
            mf_sum.append(sum_mf)
            
            if (sum_adult == 0) and (sum_juv == 0) and (sum_mf == 0):
                meta_popdict[mpop][npop] = {}
                #break

            if density_dependence: # density dependent mortality
                Km = 100  #carrying capacity
                bm = 1.0/Km #Km is carrying capacity
                Sm = 0.99 #maximum survival
                am = (Sm*bm)/math.exp(-1)
                mort_A = am * sum_adult * math.exp(-bm * sum_adult) #Ricker fx
                mort_J = am * sum_juv * math.exp(-bm * sum_juv) #Ricker fx
                mort_M = am * sum_mf * math.exp(-bm * sum_mf) #Ricker fx
                fecund = fecundity #fecundity per month
            else: # fixed mortality, actually prob of surviving
                mort_A = mortalityAdult #year
                mort_J = mortalityJuv #month
                mort_M = mortalityMF #month
                fecund = fecundity #fecundity per month

           #temp for J_12 > A_1
            #tempA1 = random.sample(meta_popdict[mpop][npop]["J_12"], int(round(len(meta_popdict[mpop][npop]["J_12"])*mort_J)))
            #for subl in tempA1:
            #    subl[0] += 1
            #first year adults
            #meta_popdict[mpop][npop]["A_1"].append(tempA1)
            # juvenilles
            meta_popdict[mpop][npop]["J_12"] = random.sample(meta_popdict[mpop][npop]["J_11"], int(round(len(meta_popdict[mpop][npop]["J_11"])*mort_J)))
            meta_popdict[mpop][npop]["J_11"] = random.sample(meta_popdict[mpop][npop]["J_10"], int(round(len(meta_popdict[mpop][npop]["J_10"])*mort_J)))
            meta_popdict[mpop][npop]["J_10"] = random.sample(meta_popdict[mpop][npop]["J_9"], int(round(len(meta_popdict[mpop][npop]["J_9"])*mort_J)))
            meta_popdict[mpop][npop]["J_9"] = random.sample(meta_popdict[mpop][npop]["J_8"], int(round(len(meta_popdict[mpop][npop]["J_8"])*mort_J)))
            meta_popdict[mpop][npop]["J_8"] = random.sample(meta_popdict[mpop][npop]["J_7"], int(round(len(meta_popdict[mpop][npop]["J_7"])*mort_J)))
            meta_popdict[mpop][npop]["J_7"] = random.sample(meta_popdict[mpop][npop]["J_6"], int(round(len(meta_popdict[mpop][npop]["J_6"])*mort_J)))
            meta_popdict[mpop][npop]["J_6"] = random.sample(meta_popdict[mpop][npop]["J_5"], int(round(len(meta_popdict[mpop][npop]["J_5"])*mort_J)))
            meta_popdict[mpop][npop]["J_5"] = random.sample(meta_popdict[mpop][npop]["J_4"], int(round(len(meta_popdict[mpop][npop]["J_4"])*mort_J)))
            meta_popdict[mpop][npop]["J_4"] = random.sample(meta_popdict[mpop][npop]["J_3"], int(round(len(meta_popdict[mpop][npop]["J_3"])*mort_J)))
            meta_popdict[mpop][npop]["J_3"] = random.sample(meta_popdict[mpop][npop]["J_2"], int(round(len(meta_popdict[mpop][npop]["J_2"])*mort_J)))
            meta_popdict[mpop][npop]["J_2"] = random.sample(meta_popdict[mpop][npop]["J_1"], int(round(len(meta_popdict[mpop][npop]["J_1"])*mort_J)))
            # microfilaria
            meta_popdict[mpop][npop]["MF_12"] = random.sample(meta_popdict[mpop][npop]["MF_11"], int(round(len(meta_popdict[mpop][npop]["MF_11"])*mort_M)))
            meta_popdict[mpop][npop]["MF_11"] = random.sample(meta_popdict[mpop][npop]["MF_10"], int(round(len(meta_popdict[mpop][npop]["MF_10"])*mort_M)))
            meta_popdict[mpop][npop]["MF_10"] = random.sample(meta_popdict[mpop][npop]["MF_9"], int(round(len(meta_popdict[mpop][npop]["MF_9"])*mort_M)))
            meta_popdict[mpop][npop]["MF_9"] = random.sample(meta_popdict[mpop][npop]["MF_8"], int(round(len(meta_popdict[mpop][npop]["MF_8"])*mort_M)))
            meta_popdict[mpop][npop]["MF_8"] = random.sample(meta_popdict[mpop][npop]["MF_7"], int(round(len(meta_popdict[mpop][npop]["MF_7"])*mort_M)))
            meta_popdict[mpop][npop]["MF_7"] = random.sample(meta_popdict[mpop][npop]["MF_6"], int(round(len(meta_popdict[mpop][npop]["MF_6"])*mort_M)))
            meta_popdict[mpop][npop]["MF_6"] = random.sample(meta_popdict[mpop][npop]["MF_5"], int(round(len(meta_popdict[mpop][npop]["MF_5"])*mort_M)))
            meta_popdict[mpop][npop]["MF_5"] = random.sample(meta_popdict[mpop][npop]["MF_4"], int(round(len(meta_popdict[mpop][npop]["MF_4"])*mort_M)))
            meta_popdict[mpop][npop]["MF_4"] = random.sample(meta_popdict[mpop][npop]["MF_3"], int(round(len(meta_popdict[mpop][npop]["MF_3"])*mort_M)))
            meta_popdict[mpop][npop]["MF_3"] = random.sample(meta_popdict[mpop][npop]["MF_2"], int(round(len(meta_popdict[mpop][npop]["MF_2"])*mort_M)))
            meta_popdict[mpop][npop]["MF_2"] = random.sample(meta_popdict[mpop][npop]["MF_1"], int(round(len(meta_popdict[mpop][npop]["MF_1"])*mort_M)))
            
            print("pass")
           # birth of MF_1
            mf1 = []
            for i in range(1, 8):
                for ad in meta_popdict[mpop][npop]["A_{}".format(i)]:
                    births = np.random.poisson(fecund)
                    n = 0
                    while n < births:
                        newmf1 = copy.copy(ad)
                        age_class = random.randint(1,8)
                        wb_parent2 = meta_popdict[mpop][npop]["A_{}".format(age_class)][random.randint(0,len(meta_popdict[mpop][npop]["A_{}".format(age_class)]))] #random adult
                        print(age_class, wb_parent2)
                        while len(wb_parent2) is 0:
                             age_class = random.randint(1,8)
                             wb_parent2 = meta_popdict[mpop][npop]["A_{}".format(age_class)][random.randint(0,len(meta_popdict[mpop][npop]["A_{}".format(age_class)]))] #random adult
                        for p in range(2, len(newmf1)):
                            newmf1[p] = newmf1[p] + wb_parent2[p]
                        mf1.append(newmf1)
                        n += 1
           #mf = sum(mf1, [])
            mf = mf1
            print(mf)
           # recombination in new mf
            num_recomb = []
            for i, bp in enumerate(basepairs):
                num_recomb.append(np.random.binomial(2*len(mf), (bp-1) * recombination[i]))
            if sum(num_recomb) != 0:
                mf = recombination_fx(mf, num_recomb, basepairs)

           #makes diploid
            for m in mf:
                for item in m:
                    if len(item) == 4:
                        item = [item[random.randint(0, 1)], item[random.randint(2, 3)]]

           #mutation in new mf
            num_muts = []
            for i in range(0, len(basepairs)):
                if recombination[i] == 0:
                    num_muts.append(np.random.binomial(len(mf), basepairs[i] * mutation[i]))
                else:
                    num_muts.append(np.random.binomial(2*len(mf), basepairs[i] * mutation[i]))
            if sum(num_muts) != 0:
                mf = mutation_fx(mf, num_muts, basepairs)
            meta_popdict[mpop][npop]["MF_1"] = mf

    return mf, meta_popdict, len(meta_popdict[mpop].keys()), sum(mf_sum)


if __name__ == '__main__':

