#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:17:00 2016
things to do: 1) add multiple loci 2) make discrete case for multiple villages
@author: scott
"""
import numpy as np
from collections import defaultdict
from math import pi, exp, sqrt
import math
import subprocess, re, random
from sklearn.metrics.pairwise import euclidean_distances
import copy
prev=[.1]
hostpopsize=[100]


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
def migration_matrix(villages):
    '''calculated a migration matrix for use with ms (hudson 2000) from euclidian
    distances. The value of 2Nm is weighted by the distance from the next population
    as an exponential random variable. The highest this can be is 2Nm/1'''  
    #distance_m is list [1000,1500] such that distance_m[0] is between 1 & 2
    #distance_m[1] is between 1 & 3 etc ...
    #villages is int 3; such that len(distance_m) == ((villages)*(villages-1)/2)
    m = 0.0001 #migration rate
    Ne = 10000 #effective size
    distance_m = [1000] #distance in meters between 2 villages
    if villages > 4: #cant figure out how to pythonically increase the villages without just if loops
        raise ValueError("only handles 4 villages ATM")
    elif villages < 4: #check that number of distances is appropriate for number of villages
        if len(distance_m) != ((villages)*(villages-1)/2): raise ValueError("there are not adequate pairwise comparisons")
        mig = [] #blank migration list
        for meters in distance_m:
            mig.append((m)/(np.random.exponential(meters)))
        if villages == 2:
            M1=4*Ne*mig[0]
            return "{} {} {} {}".format(0,M1,0,M1) #mig_matrix is symmetrical 
        elif villages == 3: 
            M1=4*Ne*mig[0]
            M2=4*Ne*mig[1]
            M3=4*Ne*mig[2]      
            return "{} {} {} {} {} {} {} {} {}".format(0,M1,M2,M1,0,M3,M2,M3,0) #mig_matrix is symmetrical
        elif villages == 4:
            M1=4*Ne*mig[0]
            M2=4*Ne*mig[1]
            M3=4*Ne*mig[2]
            M4=4*Ne*mig[3]       
            M5=4*Ne*mig[4]
            M6=4*Ne*mig[5]
            return "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(0,M1,M2,M3,M1,0,M4,M5,M2,M4,0,M6,M3,M5,M6,0)

def ms_outcall(worm_popsize):
    '''external call to ms (Hudson 2000) using the migration matrix calculated
above from the euclidian distances. This function will then read in the stdout from
ms. The continious location of Seg given by ms is transformed to discrete coordinates. 
It then builds a list of lists from the discrete coordinates which is then pointer
referenced in the dictionary recording the metapopulation information
    metapop_init = [100,100] #2 metapops
    theta = [metapop1,metapop2]    
    ''' 
    theta = [[5,2]] #4Neu expected number of mutations; [locus1],[locus2]
    rho = [0,100] #4Ner expected number of recombination of a locus per generation; [rho_1, rho_2]    
    seqL = [13000,200000]    
    theta_anc = [12,9] #theta for ancestral all pops
    ta = 0.05 #time for ancestral
    t12 = 0.05 #time for joining/splitting pop 1 into 2
    t23 = 0.05 #time for joining/splitting pop 2 into 3
    t34 = 0.05 #time for joining/splitting pop 3 into 4
    ploidy = 1    
    popsize = worm_popsize
    hap_pop = defaultdict(list) #list intialization for recording haplotypes
    
    for i in range(len(rho)):
        if rho[i] is not 0:
            ploidy = 2
        else:
            ploidy = 1
            
        popsize[:] = [x * ploidy for x in popsize]
        total_inds = sum(popsize) # how many
            
        if len(worm_popsize) == 1: #ms set up for just 1 village
            mscmd = "scrm {} 1 -t {} -r {} {} -eN {} {}".format(total_inds,theta[0][i],rho[i],seqL[i],ta,float(theta_anc[i]/theta[0][i]))        
        else: #ms setup for >1 villages
            num_subpops = len(popsize) #-I num_pops
            sub_pop = " ".join(map(str,popsize))#-I X i j ...
            if len(worm_popsize) == 2:        
                mscmd = "scrm {} 1 -t {} -r {} {} -I {} {} -n 1 {} -n 2 {} -ma {} -ej {} 1 2 -eN {} {}".format(total_inds,theta[0][i],rho[i],seqL[i],num_subpops,sub_pop,1,float(theta[1][i])/theta[0][i],migration_matrix(len(worm_popsize)),t12,ta,theta_anc[i])       
            elif len(worm_popsize)==3:
                mscmd = "scrm {} 1 -t {} -r {} {} -I {} {} -n 1 {} -n 2 {} -n 3 {} -ma {} -ej {} 1 2 -ej {} 2 3 -eN {} 1".format(total_inds,theta[0][i],rho[i],seqL[i],num_subpops,sub_pop,1,float(theta[1][i])/theta[0][i],float(theta[2][i])/theta[0][i],migration_matrix(len(worm_popsize)),t12,t23,ta,theta_anc[i])       
            elif len(worm_popsize)==4:
                mscmd = "scrm {} 1 -t {} -r {} {} -I {} {} -n 1 {} -n 2 {} -n 3 {} -n4 {} -ma {} -ej {} 1 2 -ej {} 2 3 -ej {} 3 4 -eN {} 1".format(total_inds,theta[0][i],rho[i],seqL[i],num_subpops,sub_pop,1,float(theta[1][i])/theta[0][i],float(theta[2][i])/theta[0][i],float(theta[3][i])/theta[0][i],migration_matrix(len(worm_popsize)),t12,t23,t34,ta,theta_anc[i])       
        print mscmd
        proc = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
        proc.wait()
    
        #parses ms output from stdout  
   
        for line in iter(proc.stdout.readline,''):
            if line.startswith("positions"):
                positions = map(int,line.strip().split())               
                for line in iter(proc.stdout.readline,''):
                    hap = line.strip()
                    #print hap                    
                    hap_pop["locus"+str(i)].append([m.start() for m in re.finditer("1", hap)]) #for non-recombining and infinite sites
        rel_pos = set([item for sublist in hap_pop["locus"+str(i)] for item in sublist])
    return hap_pop
    

def trans_init(prev,hostpopsize):
    '''initializes locations for above host infections within a village'''
    #set dispersal parameters
    sigma = 2  
    dispersal = 4*pi*sigma**2
    size = 1
    mu = 100
    metan = 1 #this set first village
    transmission_mat = AutoVivification() #sets dictionary     
    for i,j in zip(prev,hostpopsize):    
        x1 = np.random.negative_binomial(size,size/float((size+mu)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
        x2 = np.random.negative_binomial(size,size/float((size+mu)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
        X = np.vstack((x1,x2)).T    
        dist = euclidean_distances(X,X) #from sklearn 
        for pop in range(0,len(x1)):
            transmission_mat["meta_{}".format(metan)]["pop_{}".format(pop+1)] = ["{},{}".format(x1[pop],x2[pop]), [1 if item < dispersal else 0 for item in dist[pop,:]]]
        metan += 1
    return transmission_mat, dispersal    
  
def worm_burden(prev,hostpopsize):
    '''
    hostpopsize is list of population in each village, prev is list of percent inf per village
    mu = [5,1] #avg_burden, average number of adult female worms in infections;
    size = [1,50] #size = dispersion, size parameter for negative binomial distribution
    hap_pop is from ms_outcall    
    '''  
    # number of adult worms per infection    
    num_inf = []
    for i,j in zip(prev,hostpopsize):
        num_inf.append(round(i*j))
    pops = len(hostpopsize)
    # set parameters
    mu = np.random.uniform(1,5,pops) # average burden
    size = np.random.uniform(1,50,pops) # dispersion, how much variance in avg burden
    
    # worm burden    
    pop_init=[]
    for i,j,k in zip(mu,size,num_inf):
        wb_burden = np.random.negative_binomial(i,i/float(i+j),k) # number of successes, prob of success (size/size+mu),number of values to return
        pop_init.append(np.array(wb_burden).tolist())
    
    # call to function msoutcall to create hap_pop list
    worm_popsize = []
    for meta in pop_init:    
        worm_popsize.append(sum(meta))
    hap_pop = ms_outcall(worm_popsize)  
        
    # meta_popdict[meta1][pop][age_stage][list]
    meta_popdict = AutoVivification() #intialize dictionary  
    #initialize counters
    pop = 1
    meta = 1
    k = 0 #count lines of haplotypes and assigns them in order to populations
    kd = 0
    locus_x = 1 #autosomal loci number
    for metapop in pop_init:
        for wb_a1 in metapop: #all infections are initialized from A1. Prob need a burn in to get away from initialization
            j = 0 #initial counter       
            meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"]=[]
            while j < wb_a1:
                dipset=[np.random.uniform(), k]
                for i in range(locus_x):                    
                    dipset.append([kd,kd+1])
                #meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"].append([np.random.uniform(), k, dipset])        
                meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"].append(dipset)
                j += 1 #counts the A1 in population
                k += 1 #counts the haps in hap_pop[]
                kd += 2
            pop += 1 #advances the pop counter
        meta += 1 #advances the meta counter
        pop = 1 #resets pop counter for new meta populations aka village
    transmission_mat,dispersal=trans_init(prev,hostpopsize)
    return meta_popdict,hap_pop,transmission_mat,dispersal
    
 
def maturation(meta_popdict, hap_pop, month, prev_t):
    '''actual life cycle stage, including births/deaths of worms'''    
    hostmort = 1/70.00 # hostmortality 1 death per 70 years
    # stat collection and counters    
    rzero_freq = []
    hap_freq = []
    sum_adult_vil = [] #of adult worms in each village
    sum_juv_vil = [] #of juv worms in each village
    sum_mf_vil = [] #of MF in each village
    adult_vil = [] #of each adult stage
    juv_vil = [] #of each juv stage
    mf_vil = [] #of each mf stage
    for mpop in meta_popdict.keys():
        adult_vil.append({key: sum(len(dct[key]) for dct in meta_popdict[mpop].values()) for key in ['A_1','A_2','A_3','A_4','A_5','A_6','A_7','A_8']})
        juv_vil.append({key: sum(len(dct[key]) for dct in meta_popdict[mpop].values()) for key in ['J_1','J_2','J_3','J_4','J_5','J_6','J_7','J_8','J_9','J_10','J_11','J_12']}) 
        mf_vil.append({key: sum(len(dct[key]) for dct in meta_popdict[mpop].values()) for key in ['MF_1','MF_2','MF_3','MF_4','MF_5','MF_6','MF_7','MF_8','MF_9','MF_10','MF_11','MF_12']})
        sum_adult_vil.append(sum(adult_vil.values()))
        sum_juv_vil.append(sum(juv_vil.values()))
        sum_mf_vil.append(sum(mf_vil.values()))
        
    # Locus info, for multiple loci this should be a list
    bp = [13000,200000]
    pmut = [7.6E-8,2.9E-9]
    precomb = [0,7.6E-8]
    frac_gen = 0.125 #gen is 8 months so 1 month is .125 of a generation

    #since this month to month 
    if month%12 is 0: #this denotes 1 year has passed so adults mature to next age class 
        for mpop in meta_popdict.keys(): #villages
            r = random.uniform(0,1) 
            if r < hostmort:
                meta_popdict, prevt_t = hostcured(meta_popdict, mpop, prev_t)      
            for npop in meta_popdict[mpop].keys(): #inf
                #count individuals in each class for density dependent calculations
#                adult = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['A_1','A_2','A_3','A_4','A_5','A_6','A_7','A_8']}
#                juv = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['J_1','J_2','J_3','J_4','J_5','J_6','J_7','J_8','J_9','J_10','J_11','J_12']}
#                mf = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['MF_1','MF_2','MF_3','MF_4','MF_5','MF_6','MF_7','MF_8','MF_9','MF_10','MF_11','MF_12']}
#                sum_adult = sum(adult_vil.values())
#                sum_juv = sum(juv_vil.values())
#                sum_mf = sum(mf_vil.values())
            
            # density dependent mortality    
            #    Km = 100
            #    bm = 1.0/Km #K is carrying capacity
            #    Sm = 0.99
            #    am = (Sm*bm)/math.exp(-1) #where S is maximum survival
            #    mort_A = am * sum_adults * math.exp(-bm * sum_adult) #Ricker fx 
            #    mort_J = am * sum_juv * math.exp(-bm * sum_juv) #Ricker fx
            #    mort_M = am * sum_mf * math.exp(-bm * sum_mf) #Ricker fx
            
                # fixed mortality, actually prob of surviving  
                mort_A = .88 #year
                mort_J = .86 #month
                mort_M = .90 #month
            
            # density dep fecundity   
            #    Kf = 100
            #    bf = 1.0/Kf #K is carrying capacity
            #    Sf = 0.99
            #    af = (Sf*bf)/math.exp(-1) #where S is maximum survival
            #    fecund = 20 * af * sum_adult * math.exp(-bf * sum_adult) #Ricker fx 
                fecund = 20   #fixed fecundity
                
                # update meta_popdict
                meta_popdict[mpop][npop]["A_8"] = random.sample(meta_popdict[mpop][npop]["A_7"],int(round(len(meta_popdict[mpop][npop]["A_7"])*mort_A)))
                meta_popdict[mpop][npop]["A_7"] = random.sample(meta_popdict[mpop][npop]["A_6"],int(round(len(meta_popdict[mpop][npop]["A_6"])*mort_A)))
                meta_popdict[mpop][npop]["A_6"] = random.sample(meta_popdict[mpop][npop]["A_5"],int(round(len(meta_popdict[mpop][npop]["A_5"])*mort_A)))
                meta_popdict[mpop][npop]["A_5"] = random.sample(meta_popdict[mpop][npop]["A_4"],int(round(len(meta_popdict[mpop][npop]["A_4"])*mort_A)))
                meta_popdict[mpop][npop]["A_4"] = random.sample(meta_popdict[mpop][npop]["A_3"],int(round(len(meta_popdict[mpop][npop]["A_3"])*mort_A)))
                meta_popdict[mpop][npop]["A_3"] = random.sample(meta_popdict[mpop][npop]["A_2"],int(round(len(meta_popdict[mpop][npop]["A_2"])*mort_A)))
                meta_popdict[mpop][npop]["A_2"] = random.sample(meta_popdict[mpop][npop]["A_1"],int(round(len(meta_popdict[mpop][npop]["A_1"])*mort_A)))
                # first year adults, later add from different age classes with maturity an evolvable trait            
                meta_popdict[mpop][npop]["A_1"] = random.sample(meta_popdict[mpop][npop]["J_12"],int(round(len(meta_popdict[mpop][npop]["J_12"])*mort_J)))               
                # juvenille
                meta_popdict[mpop][npop]["J_12"] = random.sample(meta_popdict[mpop][npop]["J_11"],int(round(len(meta_popdict[mpop][npop]["J_11"])*mort_J)))
                meta_popdict[mpop][npop]["J_11"] = random.sample(meta_popdict[mpop][npop]["J_10"],int(round(len(meta_popdict[mpop][npop]["J_10"])*mort_J)))
                meta_popdict[mpop][npop]["J_10"] = random.sample(meta_popdict[mpop][npop]["J_9"],int(round(len(meta_popdict[mpop][npop]["J_9"])*mort_J)))
                meta_popdict[mpop][npop]["J_9"] = random.sample(meta_popdict[mpop][npop]["J_8"],int(round(len(meta_popdict[mpop][npop]["J_8"])*mort_J)))
                meta_popdict[mpop][npop]["J_8"] = random.sample(meta_popdict[mpop][npop]["J_7"],int(round(len(meta_popdict[mpop][npop]["J_7"])*mort_J)))
                meta_popdict[mpop][npop]["J_7"] = random.sample(meta_popdict[mpop][npop]["J_6"],int(round(len(meta_popdict[mpop][npop]["J_6"])*mort_J)))
                meta_popdict[mpop][npop]["J_6"] = random.sample(meta_popdict[mpop][npop]["J_5"],int(round(len(meta_popdict[mpop][npop]["J_5"])*mort_J)))
                meta_popdict[mpop][npop]["J_5"] = random.sample(meta_popdict[mpop][npop]["J_4"],int(round(len(meta_popdict[mpop][npop]["J_4"])*mort_J)))
                meta_popdict[mpop][npop]["J_4"] = random.sample(meta_popdict[mpop][npop]["J_3"],int(round(len(meta_popdict[mpop][npop]["J_3"])*mort_J)))
                meta_popdict[mpop][npop]["J_3"] = random.sample(meta_popdict[mpop][npop]["J_2"],int(round(len(meta_popdict[mpop][npop]["J_2"])*mort_J)))
                meta_popdict[mpop][npop]["J_2"] = random.sample(meta_popdict[mpop][npop]["J_1"],int(round(len(meta_popdict[mpop][npop]["J_1"])*mort_J)))
                # microfilaria                
                meta_popdict[mpop][npop]["MF_12"] = random.sample(meta_popdict[mpop][npop]["MF_11"],int(round(len(meta_popdict[mpop][npop]["MF_11"])*mort_M)))
                meta_popdict[mpop][npop]["MF_11"] = random.sample(meta_popdict[mpop][npop]["MF_10"],int(round(len(meta_popdict[mpop][npop]["MF_10"])*mort_M)))
                meta_popdict[mpop][npop]["MF_10"] = random.sample(meta_popdict[mpop][npop]["MF_9"],int(round(len(meta_popdict[mpop][npop]["MF_9"])*mort_M)))
                meta_popdict[mpop][npop]["MF_9"] = random.sample(meta_popdict[mpop][npop]["MF_8"],int(round(len(meta_popdict[mpop][npop]["MF_8"])*mort_M)))
                meta_popdict[mpop][npop]["MF_8"] = random.sample(meta_popdict[mpop][npop]["MF_7"],int(round(len(meta_popdict[mpop][npop]["MF_7"])*mort_M)))
                meta_popdict[mpop][npop]["MF_7"] = random.sample(meta_popdict[mpop][npop]["MF_6"],int(round(len(meta_popdict[mpop][npop]["MF_6"])*mort_M)))
                meta_popdict[mpop][npop]["MF_6"] = random.sample(meta_popdict[mpop][npop]["MF_5"],int(round(len(meta_popdict[mpop][npop]["MF_5"])*mort_M)))
                meta_popdict[mpop][npop]["MF_5"] = random.sample(meta_popdict[mpop][npop]["MF_4"],int(round(len(meta_popdict[mpop][npop]["MF_4"])*mort_M)))
                meta_popdict[mpop][npop]["MF_4"] = random.sample(meta_popdict[mpop][npop]["MF_3"],int(round(len(meta_popdict[mpop][npop]["MF_3"])*mort_M)))
                meta_popdict[mpop][npop]["MF_3"] = random.sample(meta_popdict[mpop][npop]["MF_2"],int(round(len(meta_popdict[mpop][npop]["MF_2"])*mort_M)))
                meta_popdict[mpop][npop]["MF_2"] = random.sample(meta_popdict[mpop][npop]["MF_1"],int(round(len(meta_popdict[mpop][npop]["MF_1"])*mort_M)))

                #add generation to A_1 since they are new from Juv12 and 1 year has passed so all are moving to A2
                for subl in meta_popdict[mpop][npop]["A_1"]:
                    subl[0] += 1

                # birth of MF_1
                mf1 = []                
                for i in range(1,8):
                    for j in meta_popdict[mpop][npop]["A_{}".format(i)]:
                        births = np.random.poisson(fecund)                        
                        while n < births:                         
                            newmf1 = copy.copy(j) #j[2][i] is vector of diploid
                            wb_parent2 = meta_popdict[mpop][npop]["A_{}".format(random.randint(1,8))] #                           
                            for p in range(2,len(newmf1)):                             
                                newmf1[p] = newmf1[p] + wb_parent2[p]
                            mf1.append(newmf1)
                            n += 1
                ???? mf = sum(mf1,[])
                
                # recombination in new mf
                num_recomb = []
                for i in range(len(bp)):
                    num_recomb.append(np.random.binomial(2*len(mf), bp[i] * precomb[i])) 
                if sum(num_recomb) != 0:
                    mf,hap_pop = recombination(mf,hap_pop,num_recomb)
                    
                #makes diploid
                for m in mf:
                    for item in m:
                        if len(item) == 4:
                            item = [item[random.randint(0,1)],item[random.randint(2,3)]] 
                            
                #mutation in new mf  
                num_muts = []
                for i in range(0,len(bp)):
                    if precomb[i] == 0:
                        num_muts.append(np.random.binomial(len(mf), bp[i] * pmut[i]))
                    else:
                        num_muts.append(np.random.binomial(2*len(mf), bp[i] * pmut[i]))               
                if sum(num_muts) != 0:
                    mf,hap_pop = mutation(mf,hap_pop,num_muts)
                meta_popdict[mpop][npop]["MF_1"] = mf
                
                #hap_freq.append([i[1]for i in mf])
                
    else: #a year has not passed on months, juveniles and MF move to next age class
        for mpop in meta_popdict.keys():        
            for npop in meta_popdict[mpop].keys(): #inf                
                #count individuals in each class for density dependent calculations
#                adult = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['A_1','A_2','A_3','A_4','A_5','A_6','A_7','A_8']}
#                juv = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['J_1','J_2','J_3','J_4','J_5','J_6','J_7','J_8','J_9','J_10','J_11','J_12']}
#                mf = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['MF_1','MF_2','MF_3','MF_4','MF_5','MF_6','MF_7','MF_8','MF_9','MF_10','MF_11','MF_12']}
#                sum_adult = sum(adult_vil.values())
#                sum_juv = sum(juv_vil.values())
#                sum_mf = sum(mf_vil.values())
            
            # density dependent mortality    
            #    Km = 100
            #    bm = 1.0/Km #K is carrying capacity
            #    Sm = 0.99
            #    am = (Sm*bm)/math.exp(-1) #where S is maximum survival
            #    mort_A = am * sum_adults * math.exp(-bm * sum_adult) #Ricker fx 
            #    mort_J = am * sum_juv * math.exp(-bm * sum_juv) #Ricker fx
            #    mort_M = am * sum_mf * math.exp(-bm * sum_mf) #Ricker fx
            
                # fixed mortality, actually prob of surviving  
                mort_A = .88 #year
                mort_J = .86 #month
                mort_M = .90 #month
            
            # density dep fecundity   
            #    Kf = 100
            #    bf = 1.0/Kf #K is carrying capacity
            #    Sf = 0.99
            #    af = (Sf*bf)/math.exp(-1) #where S is maximum survival
            #    fecund = 20 * af * sum_adult * math.exp(-bf * sum_adult) #Ricker fx 
                fecund = 20   #fixed fecundity

                # juvenilles                
                meta_popdict[mpop][npop]["J_12"].append(random.sample(meta_popdict[mpop][npop]["J_11"],int(round(len(meta_popdict[mpop][npop]["J_11"])*mort_J))))
                meta_popdict[mpop][npop]["J_11"] = random.sample(meta_popdict[mpop][npop]["J_10"],int(round(len(meta_popdict[mpop][npop]["J_10"])*mort_J)))
                meta_popdict[mpop][npop]["J_10"] = random.sample(meta_popdict[mpop][npop]["J_9"],int(round(len(meta_popdict[mpop][npop]["J_9"])*mort_J)))
                meta_popdict[mpop][npop]["J_9"] = random.sample(meta_popdict[mpop][npop]["J_8"],int(round(len(meta_popdict[mpop][npop]["J_8"])*mort_J)))
                meta_popdict[mpop][npop]["J_8"] = random.sample(meta_popdict[mpop][npop]["J_7"],int(round(len(meta_popdict[mpop][npop]["J_7"])*mort_J)))
                meta_popdict[mpop][npop]["J_7"] = random.sample(meta_popdict[mpop][npop]["J_6"],int(round(len(meta_popdict[mpop][npop]["J_6"])*mort_J)))
                meta_popdict[mpop][npop]["J_6"] = random.sample(meta_popdict[mpop][npop]["J_5"],int(round(len(meta_popdict[mpop][npop]["J_5"])*mort_J)))
                meta_popdict[mpop][npop]["J_5"] = random.sample(meta_popdict[mpop][npop]["J_4"],int(round(len(meta_popdict[mpop][npop]["J_4"])*mort_J)))
                meta_popdict[mpop][npop]["J_4"] = random.sample(meta_popdict[mpop][npop]["J_3"],int(round(len(meta_popdict[mpop][npop]["J_3"])*mort_J)))
                meta_popdict[mpop][npop]["J_3"] = random.sample(meta_popdict[mpop][npop]["J_2"],int(round(len(meta_popdict[mpop][npop]["J_2"])*mort_J)))
                meta_popdict[mpop][npop]["J_2"] = random.sample(meta_popdict[mpop][npop]["J_1"],int(round(len(meta_popdict[mpop][npop]["J_1"])*mort_J)))
                # microfilaria                
                meta_popdict[mpop][npop]["MF_12"] = random.sample(meta_popdict[mpop][npop]["MF_11"],int(round(len(meta_popdict[mpop][npop]["MF_11"])*mort_M)))
                meta_popdict[mpop][npop]["MF_11"] = random.sample(meta_popdict[mpop][npop]["MF_10"],int(round(len(meta_popdict[mpop][npop]["MF_10"])*mort_M)))
                meta_popdict[mpop][npop]["MF_10"] = random.sample(meta_popdict[mpop][npop]["MF_9"],int(round(len(meta_popdict[mpop][npop]["MF_9"])*mort_M)))
                meta_popdict[mpop][npop]["MF_9"] = random.sample(meta_popdict[mpop][npop]["MF_8"],int(round(len(meta_popdict[mpop][npop]["MF_8"])*mort_M)))
                meta_popdict[mpop][npop]["MF_8"] = random.sample(meta_popdict[mpop][npop]["MF_7"],int(round(len(meta_popdict[mpop][npop]["MF_7"])*mort_M)))
                meta_popdict[mpop][npop]["MF_7"] = random.sample(meta_popdict[mpop][npop]["MF_6"],int(round(len(meta_popdict[mpop][npop]["MF_6"])*mort_M)))
                meta_popdict[mpop][npop]["MF_6"] = random.sample(meta_popdict[mpop][npop]["MF_5"],int(round(len(meta_popdict[mpop][npop]["MF_5"])*mort_M)))
                meta_popdict[mpop][npop]["MF_5"] = random.sample(meta_popdict[mpop][npop]["MF_4"],int(round(len(meta_popdict[mpop][npop]["MF_4"])*mort_M)))
                meta_popdict[mpop][npop]["MF_4"] = random.sample(meta_popdict[mpop][npop]["MF_3"],int(round(len(meta_popdict[mpop][npop]["MF_3"])*mort_M)))
                meta_popdict[mpop][npop]["MF_3"] = random.sample(meta_popdict[mpop][npop]["MF_2"],int(round(len(meta_popdict[mpop][npop]["MF_2"])*mort_M)))
                meta_popdict[mpop][npop]["MF_2"] = random.sample(meta_popdict[mpop][npop]["MF_1"],int(round(len(meta_popdict[mpop][npop]["MF_1"])*mort_M)))

                # birth of MF_1
      
                # recombination in new mf

                
                #makes diploid


                
                # mutation in new MF_1 population         


                #hap_freq.append([i[1]for i in mf])

    #this will calculate for the enitre meta and each inf
    #rzero_freq = sum(rzero_freq,[])
    #hap_freq = sum(hap_freq,[])
    
    return meta_popdict, hap_pop, [sum_adult_vil,sum_juv_vil,sum_mf_vil],{x:rzero_freq.count(x) for x in rzero_freq},{y:hap_freq.count(y) for y in hap_freq}     

def mutation(mf,hap_pop,num_muts,bp):
   '''this is run every time the prob of a mutation is true'''      
   #flatten list to check for colisions: set([item for sublist in hap_pop["locus1"] for item in sublist])   
   for locus in range(len(num_muts)): #number of muts
       muts = 0       
       while muts < locus: #keep going until all mutations are  assigned
           mut_mf = random.randrange(0,len(mf)) #choose random index in mf list
               #add new sequence to hap_pop
           new_hap = copy.copy(hap_pop["locus" + str(locus)][mf[mut_mf][locus+1]])
           #random position            
           new_allele = random.randint(0,bp[locus])
           new_hap.append(new_allele)
           hap_pop["locus" + str(locus)].append(new_hap.sort())
           if len(mf[mut_mf][locus+1]) > 1:
               mf[mut_mf][locus+1][random.randint(0,1)] = len(hap_pop["locus" + str(locus)])-1 #last entry is new_hap
           else:
               mf[mut_mf][locus+1]= len(hap_pop["locus" + str(locus)])-1
           muts += 1
   return mf, hap_pop

def recombination(mf,hap_pop,num_recomb):
       '''this is run every time the prob of a recombination is true'''      
   for locus in range(len(num_recomb)): #number of muts
       recomb = 0       
       while recomb < locus: #keep going until all recombinations are  assigned
           rec_mf = random.randrange(0,len(mf)) #choose random index in mf list
               #add new sequence to hap_pop
           new_recomb = random.randint(0,3)
           if new_recomb < 2:
               #first parent; create new hap
               hap1 = mf[rec_mf][locus+1][0]
               hap2 = mf[rec_mf][locus+1][1]
               new_hap = hap1 + hap2
               hap_pop["locus" + str(locus)].append(new_hap)
               mf[rec_mf][locus+1][random.randint[0,1]] = len(hap_pop["locus" + str(locus)])-1
           else:
               #second parent
               hap1 = mf[rec_mf][locus+1][0]
               hap2 = mf[rec_mf][locus+1][1]
               new_hap = hap1 + hap2
               mf[rec_mf][locus+1][random.randint[2,3]] = len(hap_pop["locus" + str(locus)])-1
               
           recomb += 1                      
           
   return mf, hap_pop
   
def transmission(transmission_mat,mpop,meta_popdict,dispersal):    
    '''this has a continious prob using a gillepsie algorithm to get the waiting time between transmission events'''    
    #this executes if tranmission is true     
    #pick a random donating pop:
    dpop = random.choice(meta_popdict[mpop].keys())    
    while [len(meta_popdict[mpop][dpop][i]) is 0 for i in meta_popdict[mpop][dpop] if "MF" in i].count(True) is 12:
        dpop = random.choice(meta_popdict[mpop].keys())
    r = random.uniform(0,1)    
    if len(transmission_mat[mpop][dpop][1]) >= hostpopsize:
        prob_newinf = 0
    else:
        prob_newinf = (float(1)/len(transmission_mat[mpop][dpop][1]))
    #new infection
    if r < prob_newinf:
        newpop = "pop_{}".format(len(meta_popdict[mpop].keys())+1)
        meta_popdict[mpop][newpop]["J_1"] = []
        rmf = int(round(random.uniform(.51,12)))
        while len(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) is 0:
            rmf = int(round(random.uniform(.51,12)))
        dmf = random.choice(meta_popdict[mpop][dpop]["MF_{}".format(rmf)])      
        meta_popdict[mpop][newpop]["J_1"].append(dmf)
        print "new"
        transmission_mat = new_infection(transmission_mat, mpop,dpop,newpop,dispersal) #find new transmission position
    #reinfection
    else:     
        rpop = random.choice([i for i, x in enumerate(transmission_mat[mpop][dpop][1]) if x == 1]) + 1 #choose a random value that is 1 as the receiving pop        
        if meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"] is not list:
            meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"] = []
        rmf = int(round(random.uniform(.51,12)))
        while len(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) is 0:
            rmf = int(round(random.uniform(.51,12)))
        dmf = random.choice(meta_popdict[mpop][dpop]["MF_{}".format(rmf)])      
        meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"].append(dmf)
    return meta_popdict, transmission_mat

def new_infection(transmission_mat, mpop,dpop,newpop,dispersal):
    '''this is run every time a new individual is infected to rebuild the transmission prob matrix, updates'''        
    x1n = np.random.uniform(-1*sqrt(dispersal), sqrt(dispersal)) 
    x2n = np.random.uniform(-1*sqrt(dispersal), sqrt(dispersal))    
    transmission_mat[mpop][newpop] = ["{},{}".format(int(transmission_mat[mpop][dpop][0].split(",")[0]) + x1n, int(transmission_mat[mpop][dpop][0].split(",")[1]) + x2n),[]]
    for tpop in range(1,len(transmission_mat[mpop].keys())):
        dist = sqrt((float(transmission_mat[mpop]["pop_{}".format(tpop)][0].split(",")[0]) - x1n)**2 + (float(transmission_mat[mpop]["pop_{}".format(tpop)][0].split(",")[1]) - x2n)**2)        
        if dist < dispersal:
            dist = 1
        else: 
            dist = 0
        transmission_mat[mpop][newpop][1].append(dist)
        transmission_mat[mpop]["pop_{}".format(tpop)][1].append(dist)
    transmission_mat[mpop][newpop][1].append(1)   
    #if pop dies it is N in the transmission_mat, it is removed from meta_popdict
     
    return transmission_mat

def hostcured(meta_popdict, mpop, prev_t):
    dpop = random.choice(meta_popdict[mpop].keys()) #draw random host infection
    meta_popdict[mpop][dpop] = {} #blanks the dictionary
    prev_t[int(dpop.split("_")[1])] = 0 #update prev_t
    return meta_popdict, prev_t
    
def wb_sims(numberGens):
    '''this will call other functions to intialize, then functions from the life cycle
    num_gen: how long to run
    prevalance: .80 is 80%
    hostpop_size: number of potentially infected hosts'''
    f = open("sims.out",'w')
    f.write("time\tvillage\tnpops\tnumA1\tnumA2\tnumA2\tnumA4\tnumA5\tnumA6\tnumA7\tnumA8\tnumMF\ttrans_event\tRzero\tnum_uniqHapsAdult\tnum_uniqHapsMF\n")
    #initialize
    meta_popdict, hap_pop, transmission_mat, dispersal = worm_burden(prev, hostpopsize)  
    #set counters
    time_month = 0
    #mf = 100
    mpop = 1
    #set prev counter
    prev_t = [0] * hostpopsize[0]
    X = prev[0] * hostpopsize[0]
    prev_t[0:X] = [1] * X 
    #how long to run the simulation
    sim_time = numberGens * 10 #gens per year are 1.2, 12/1.2 is 10 months for 1 generation            
    while time_month <= sim_time:
        #transmission
        totalbitesvillage = 20 * 6 * 30 * hostpopsize[0] #bites per person per hour * 6 hours (atnight) * 30 days * hostpopsize
        infbites=np.random.binomial(totalbitesvillage,(prev_t*0.37))
        L3trans=np.random.binomial(infbites,(0.13248))
#        mf = (((wb_pop[2])/hostpopsize*prev_t)/235)/50 #total MF in the village
#        L3trans = infbites * (4.395*(1-math.exp((.055*(mf))/4.395))**2) * 0.13248
        for i in range(L3trans):
            meta_popdict, transmission_mat = transmission(transmission_mat,mpop,meta_popdict,dispersal)  
        #maturation
        meta_popdict, hap_pop, wb_pop, rzero_freq, hap_freq = maturation(meta_popdict,hap_pop,time_month) #update matrices    
        
        #update prev_t
        prev_t = [i>0 for i in list(len(lst) for pop in meta_popdict["meta_1"].values() for lst in pop.values())].count(True) / float(hostpopsize[0])
        #len(meta_popdict[mpop].values()) #number of pops if all are not 0 then this is prevs
        
        #sum data for month and print to out file
        f.write("%i\t%s\t%f\t%i\t%i\t%i\f%i\t%i\n" %(time_month,"meta_1",prev,wb_pop[0],wb_pop[1],wb_pop[2],np.mean(rzero_freq.values()),len(hap_freq.keys())))
        time_month += 1
            
    f.close()








