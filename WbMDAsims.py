#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:17:00 2016
requires: numpy, sklearn, scrm, R 
Anaconda will install all dependencies as well as scrm 
(e.g., conda install -c bioconda scrm)
scrm: Paul R. Staab et al ..., Bioinformatics 2015 
https://scrm.github.io/

##PROGRAM STRUCTURE
wb_sims(numberGens, prev, hostpopsize, bites_person, hours2bite)
  worm_burden(prev, hostpopsize, muWormBurden, sizeWormBurden)
    ms_outcall(*worm_popsize,villages, theta, recombination, basepairs, mutation, time2Ancestral, thetaAncestral, thetaRegional, time_join12, time_join23, time_join34)
      migration_matrix(villages, initial_migration, initial_distance_m, theta[i], basepairs[i], mutation[i])
    trans_init(prev, hostpopsize, muTrans, sigma, sizeTrans)
  maturation(meta_popdict, month, prev_t, density_dependence, mortalityHost, mortalityAdult, mortalityJuv, mortalityMF, fecundity, basepairs, mutation, recombination)
    hostdeath(meta_popdict, *mpop, prev_t)
    recombination(*mf, *num_recomb, basepairs)
    mutation(*mf, *num_mut, basepairs)
  transmission(transmission_mat, *mpop, meta_popdict, dispersal)
    new_infection(transmission_mat, mpop, dpop, newpop, dispersal)
##
@author:stsmall
"""
import numpy as np
from collections import defaultdict
from math import pi, exp, sqrt, log
import math
import subprocess, re, random
from sklearn.metrics.pairwise import euclidean_distances
import copy
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    ## initialize
    #migration_matrix
    parser.add_argument('-im', '--initial_migration',type=float,default=.0001,help="migration rate between villages/metapopulations in the model.This is strictly for initial conditions and can be changed for the forward-in-time portion")
    parser.add_argument('-idm','--initial_distance_m', type=list,help="initial_distance_m is list [1000] such that distance_m[0] is between 1 & 2")
    parser.add_argument('-v','--villages', type=int, default=1, help="sets the intial number of villages.")
    parser.add_argument('-t', '--theta',type=list,required=True,help="theta value of worm populations in 1st village, should be for the entire length of the locus")
    parser.add_argument('-bp','--basepairs',type=list,default=13000,help="length in basepairs of each locus")
    parser.add_argument('-u','--mutation',type=list,default=7.6E-8,help="expected as list, mutation rate per bp per generation for each locus")    
    #ms_outcall
    parser.add_argument('-t12', '--time_join12',type=int,default=240, help="generations until time of joining for pop1 and pop2")
    parser.add_argument('-t23', '--time_join23',type=int,default=240,help="generations until time of joining for pop2 and pop3")
    parser.add_argument('-t34', '--time_join34',type=int,default=240,help="generations until time of joining for pop3 and pop4")
    parser.add_argument('-t2a','--time2Ancestral',type=int, default=1800,help="generations until ancestral population for PNG is 1800 generations for Africa/Haiti 500 generations")
    parser.add_argument('-at','--thetaAncestral',type=int,default=344, help="ancestral theta before")
    parser.add_argument('-ar','--thetaRegional',type=int,default=23, help="regional theta")
    parser.add_argument('-r','--recombination',type=list, default=0, help="recombination for each locus, if 0 assumes haploid is 1")    
    #trans_init
    parser.add_argument('-mt', '--muTrans',type=int,default=100, help="mu for neg bino in transmission, distances between hosts")
    parser.add_argument('-st', '--sizeTrans',type=int,default=1, help="size for neg bino in transmission, distance between hosts")
    parser.add_argument('-dp', '--sigma',type=int,default=2, help="sigma, for dispersal")
    #worm_burden
    parser.add_argument('-prev', '--prev', type=list, required=True, help="prevalance of Wb in host populations; should be a list")
    parser.add_argument('-host', '--hostpopsize',type=list, required=True,help="size of host population")
    parser.add_argument('-mw', '--muWormBurden',type=list, help="mu for negative binomial in worm burden per host, add value or defined by uniform")
    parser.add_argument('-sw', '--sizeWormBurden',type=list, help="size for negative binomial in worm burden per host")
    ## wb_sims
    parser.add_argument('-ng','--numberGens', type=int, required=True, help="total number of generation to run the simulation")
    parser.add_argument('-bpp','--bites_person', type=int, default = 20, help="number of bites recieved per person per hour")
    parser.add_argument('-h2b','--hours2bite', type=int, default = 6, help="number of hours exposed to mosquitoes")
    #wblifecycles
    parser.add_argument('-ma','--mortalityAdult',type=float,default=0.8863848717161292,help="adult worms dies at rate 0.01 per month or survive at .99")
    parser.add_argument('-mj','--mortalityJuv',type=float,default=0.865880173963825,help="prob that juv makes it to adult is 0.2 or 0.8177651 per month")
    parser.add_argument('-mf','--mortalityMF',type=float,default=0.90,help="MF die at rate 0.10 per month or survive at 0.90")
    parser.add_argument('-mh','--mortalityHost',type=float,default=0.14286,help="Host death per year")
    parser.add_argument('-f','--fecundity',type=int,default=20,help="mean number of MF born to each Adult per month")
    parser.add_argument('-Dd','--density_dependence',action="store_true",help="use density dependence")
    #if args.density_dependence:
    parser.add_argument('-gtime','--generations',type=int,default=0.125,help="generation time in years")
    parser.add_argument('-hostmig','--host_migration_rates', help="list of host migration rates between villages per month")
    args = parser.parse_args()
    return args

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
 
#################################################
##Initialize genetic,hosts,transmission
#################################################
            
def migration_matrix(villages,initial_migration,initial_distance_m,theta,basepairs,mutation):
    '''creates a string that represents a migration matrix between villages (metapopulations)
    in the format determined by ms (hudson 2000).Uses euclidian distances.The value of 2Nm 
    is weighted by the distance from the next population as an exponential random variable. 
    The highest this can be is 2Nm/1
    villages:(int) number of villages/metapopulations
    initial_migration:(float) migration for ms/scrm as the fraction of subpopulation i which is made up of migrants from subpopulation j each generation
    initial_distance_m:(list,float) distance between villages in meters
    '''
    Ne = 4 * basepairs * theta / mutation
    if villages > 4: #cant figure out how to pythonically increase the villages without just using loops
        raise ValueError("only handles 4 villages ATM")
    elif villages < 4: 
        if len(initial_distance_m) != ((villages)*(villages-1)/2): #check that number of distances is appropriate for number of villages 
            raise ValueError("there are not adequate pairwise comparisons in distance_m to match villages")
        mig = [] # initiate blank migration list
        for meters in initial_distance_m:
            mig.append((initial_migration)/(np.random.exponential(meters)))
        if villages == 2:
            M1=4*Ne*mig[0] #4Nm
            return "{}".format(M1) #mig_matrix is symmetrical and island 
        elif villages == 3: 
            M1 = 4*Ne*mig[0]
            M2 = 4*Ne*mig[1]
            M3 = 4*Ne*mig[2]      
            return "{} {} {} {} {} {} {} {} {}".format(0,M1,M2,M1,0,M3,M2,M3,0) #mig_matrix is symmetrical
        elif villages == 4:
            M1 = 4*Ne*mig[0]
            M2 = 4*Ne*mig[1]
            M3 = 4*Ne*mig[2]
            M4 = 4*Ne*mig[3]       
            M5 = 4*Ne*mig[4]
            M6 = 4*Ne*mig[5]
            return "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(0,M1,M2,M3,M1,0,M4,M5,M2,M4,0,M6,M3,M5,M6,0)

def ms_outcall(worm_popsize,villages,initial_migration,initial_distance_m,theta,basepairs,mutation,recombination,time2Ancestral,thetaAncestral,thetaRegional,time_join12,time_join23,time_join34):
    '''external call to ms (Hudson 2000) or scrm. Calls the function migration_matrix() if using more than 1 village.
    This function will then read in the stdout from ms/scrm.The location of segsites are assigned a random number then scaled
    by the length of the desired sequence in basepairs. 
    worm_popsize:(list,int) how many worms to simulate
    villages:(int) number of villages/metapopulations
    theta:(list,float) list of theta values for each locus
    basepairs:(list,int) list of basepairs (lengths) for each locus
    mutation:(list,float) mutation rate for each locus as probability per base per generation
    recombination:(list,float) recombination rate for each locus as probability per base per generation 
    thetaAncestral:(float) theta for ancestral pops to the ratio of N0; e.g., 
    thetaRegional: (float) theta for regional pops to the ration of N0 e.g.,23 times larger
    time2Ancestral:(int) time in generations to the ancestral population
    time_join12:(int) time in generations for joining/splitting pop 1 into 2
    time_join23:(int) time in generations for joining/splitting pop 2 into 3
    time_join34:(int) time in generations for joining/splitting pop 3 into 4
    '''
    num_loci = len(theta) #counts number of loci
    thetaN0 = [t[0] for t in theta] #theta for first population, all other thetas are scaled from this value
    rho = ([(b-1)*r*t/m for b,r,t,m in zip(basepairs,recombination, thetaN0, mutation)]) #population recombination rate
    tA = ([b*(a/(t/m)) for b,a,t,m in zip(basepairs,time2Ancestral, thetaN0, mutation)]) #time to the ancestral population
    t12 = ([b*(a/(t/m)) for b,a,t,m in zip(basepairs,time_join12, thetaN0, mutation)]) # time to join 1 & 2
    t23 = ([b*(a/(t/m)) for b,a,t,m in zip(basepairs,time_join23, thetaN0, mutation)]) # time to join 2 & 3
    t34 = ([b*(a/(t/m)) for b,a,t,m in zip(basepairs,time_join34, thetaN0, mutation)])  # time to join 3 & 4   
    hap_pop = defaultdict(list) #list intialization for recording haplotypes

    for i in range(0,num_loci):
        if rho[i] is 0:
            ploidy = 1
        else:
            ploidy = 2
        worm_popsize[:] = [x * ploidy for x in worm_popsize]
        total_inds = sum(worm_popsize) # how many
        #positions = '%.2E' %(Decimal(str(basepairs[i])))
        
        #create ms/scrm command    
        if villages == 1: #set up for just 1 village, doesnot call migration_matrix
            mscmd = "scrm {} 1 -t {} -r {} {} -G {} -eG {} 0.0 -eN {} {} -SC abs -p {}".format(total_inds,thetaN0[i],rho[i],basepairs[i]-1,(-1/tA[i])*math.log(thetaRegional),tA[i],tA[i],thetaAncestral,12)        
        else: #ms setup for >1 villages
            num_subpops = len(worm_popsize) #-I num_pops
            sub_pop = " ".join(map(str,worm_popsize))#-I X i j ...
            if villages == 2:        
                mscmd = "scrm {} 1 -t {} -r {} {} -I {} {} {} -n 1 {} -n 2 {} -ej {} 1 2 -G {} -eG {} 0.0 -eN {} {} -SC abs -p {}".format(total_inds,theta[i],rho[i],basepairs[i]-1,num_subpops,sub_pop,migration_matrix(villages,initial_migration,initial_distance_m,thetaN0[i],basepairs[i],mutation[i]),1,float(thetaN0[i])/theta[i][1],t12[i],(-1/tA[i])*math.log(thetaRegional),tA[i],tA[i],thetaAncestral,12)       
            elif villages == 3:
                mscmd = "scrm {} 1 -t {} -r {} {} -I {} {} -n 1 {} -n 2 {} -n 3 {} -ma {} -ej {} 1 2 -ej {} 2 3 -G {} -eG {} 0.0 -eN {} {} -SC abs -p {}".format(total_inds,theta[i],rho[i],basepairs[i]-1,num_subpops,sub_pop,migration_matrix(villages,initial_migration,initial_distance_m,thetaN0[i],basepairs[i],mutation[i]),1,float(thetaN0[i])/theta[i][1],t12[i],t23[i],(-1/tA[i])*math.log(thetaRegional),tA[i],tA[i],thetaAncestral,12)       
            elif villages == 4:
                mscmd = "scrm {} 1 -t {} -r {} {} -I {} {} -n 1 {} -n 2 {} -n 3 {} -n4 {} -ma {} -ej {} 1 2 -ej {} 2 3 -ej {} 3 4 -G {} -eG {} 0.0 -eN {} {} -SC abs -p {}".format(total_inds,theta[i],rho[i],basepairs[i]-1,num_subpops,sub_pop,migration_matrix(villages,initial_migration,initial_distance_m,thetaN0[i],basepairs[i],mutation[i]),1,float(thetaN0[i])/theta[i][1],t12[i],t23[i],t34[i],(-1/tA[i])*math.log(thetaRegional),tA[i],tA[i],thetaAncestral,12)       

        #print mscmd
        proc = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
        proc.wait()
    
        #parses ms output from stdout  
        for line in iter(proc.stdout.readline,''):
            if line.startswith("positions"):
                positions = [round(i) for i in map(int,line.strip().split()[1:])]               
                for line in iter(proc.stdout.readline,''):
                    hap = line.strip()
                    hap_pop["locus"+str(i)].append(filter(lambda a: a != 0, [i*j for i,j in zip([m.start() for m in re.finditer("1", hap)],positions)])) 
    return hap_pop
    

def trans_init(prev,hostpopsize,muTrans,sigma,sizeTrans):
    '''Creates a transmission matrix for locations of infected hosts
    muTrans:(float) parameter of neg binomial
    sizeTrans:(float) parameter of neg binomial
    sigma:(float) parameter for calculating dispersal distance
    '''
    #set dispersal parameters
    dispersal = 4*pi*sigma**2
    meta_n = 1 #this set first village
    transmission_mat = AutoVivification() #creates dictionary   ### this might not be necessary if we leave as matrix type
     
    for i,j in zip(prev,hostpopsize):    
        x1 = np.random.negative_binomial(sizeTrans,sizeTrans/float((sizeTrans+muTrans)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
        x2 = np.random.negative_binomial(sizeTrans,sizeTrans/float((sizeTrans+muTrans)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
        X = np.vstack((x1,x2)).T #transposes
        dist = euclidean_distances(X,X) #from sklearn to calculate euclidian distances 
        for pop in range(0,round(hostpopsize*prev)):
            transmission_mat["meta_{}".format(meta_n)]["pop_{}".format(pop+1)] = ["{},{}".format(x1[pop],x2[pop]), [1 if item < dispersal else 0 for item in dist[pop,:]]]
            ### is it potential for rather than 0 or 1 given success of transmission that it could be exponetial or similar long-tailed distribution
            ### however that would mean change in how transmission occurs since it would not always be a success
        meta_n += 1
    return transmission_mat, dispersal    
  
def worm_burden(prev,hostpopsize,muWormBurden,sizeWormBurden,villages,initial_migration,initial_distance_m,theta,basepairs,mutation,recombination,time2Ancestral,thetaAncestral,thetaRegional,time_join12,time_join23,time_join34,muTrans,sizeTrans,sigma):
    '''
    worm_burden calls all the previous functions to create meta_popdict. 
    meta_popdict has the age/stage structure as well as the haplotypes for each locus.
    transmission_mat gives the corresponding locations of the infected hosts
    prev: (list, float) percent of hosts infected
    hostpopsize:(list, int) number of possible hosts (human population)
    muWormBurden: (list,int) avg_burden, average number of adult female worms in infections;
    sizeWormBurden:(list,int) dispersion, size parameter for negative binomial distribution   
    '''
    # number of adult worms per infection    
    num_inf = []
    for i,j in zip(prev,hostpopsize):
        num_inf.append(round(i*j))
        
    # worm burden (number of adult worms) per host    
    pop_init=[]
    for mu,size,numworms in zip(muWormBurden,sizeWormBurden,num_inf):
        wb_burden = np.random.negative_binomial(mu,mu/float(mu+size),numworms) # number of successes, prob of success (size/size+mu),number of values to return
        pop_init.append(np.array(wb_burden).tolist())
    
    # call to function msoutcall to create hap_pop list
    worm_popsize = []
    for meta in pop_init:    
        worm_popsize.append(sum(meta))
    hap_pop = ms_outcall(worm_popsize,villages,initial_migration,initial_distance_m,theta,basepairs,mutation,recombination,time2Ancestral,thetaAncestral,thetaRegional,time_join12,time_join23,time_join34)  
        
    # meta_popdict[meta1][pop][age_stage][list]
    meta_popdict = AutoVivification() #intialize dictionary  
    #initialize counters
    pop = 1 #pop counter, also infected hosts
    meta_n = 1 #village counter
    k = 0 #count lines of haploid haplotypes and assigns them in order to individuals
    kd = 0 #count lines of diploid haplotypes and assigns them in order to individuals
    locus = len(recombination) #total loci to add
    for vill in pop_init:
        for hostpop in vill: #all infections are initialized from A1. Prob need a burn in to get away from initialization
            j = 0 #tracks number of worms to add to each hostpop       
            meta_popdict["meta_" + str(meta_n)]["pop_" + str(pop)]["A_1"]=[]  #only A_1 worms
            while j < hostpop:
                haplotypes=[np.random.uniform(),[hap_pop["locus0"][k]]] #initializes info for each worm [rand_num,[],[,],[,],[,],[,]]
                for i in range(1,locus):   
                    haplotypes.append(hap_pop["locus"+str(i)][kd],hap_pop["locus"+str(i)][kd+1])
                meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"].append(haplotypes)
                j += 1 #counts the A1 in pop_init
                k += 1 #counts the haps in hap_pop[]
                kd += 2 #counts the diploid haps in hap_pop[]
            pop += 1 #advances the pop counter
        meta_n += 1 #advances the meta counter
        pop = 1 #resets pop counter for new meta populations aka village
    transmission_mat,dispersal=trans_init(prev,hostpopsize,muTrans,sigma,sizeTrans)
    return meta_popdict,transmission_mat,dispersal

############################################
##Begin simulation functions
############################################    
   
def maturation(mpop,meta_popdict,time_month,infhost,mf_mpopavg,density_dependence,mortalityHost,mortalityAdult,mortalityJuv,mortalityMF,fecundity,basepairs,mutation,recombination):
    '''life cycle stage, including births/deaths of worms
    meta_popdict: (dict) perform functions on dictionary containing vill,host,worm and genotypes
    time_month: (int) what month is it 1 ... 12
    infhost: (int) hostpopsize * prev = number of infected hosts 
    density_dependence: (Boolean) if True then use density dependencec functions to calculate mortality/fecundity
    mortalityHost:(float) prob of host dying per year, default 1/70
    mortalityAdult:(float) prob of adult surviving per year
    mortalityJuv:(float) prob of Juv surviving per month
    mortalityMF:(float) prob of MF surviving per month
    fecundity:(int) average number of MF produced per adult per month
    '''
    mpop = "meta_" + str(mpop+1)
    mf_sum = []    
    infhost_t = 0
    #since this month to month 
    if time_month%12 is 0: #this denotes 1 year has passed so adults mature to next age class 
       for npop in meta_popdict[mpop].keys(): #infected hosts as pops within villages
           if random.uniform(0,1) < mortalityHost:
               meta_popdict = hostdeath(meta_popdict, mpop)      
           
           #count individuals in each class for density dependent calculations
           adult_dd = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['A_1','A_2','A_3','A_4','A_5','A_6','A_7','A_8']}
           juv_dd = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['J_1','J_2','J_3','J_4','J_5','J_6','J_7','J_8','J_9','J_10','J_11','J_12']}
           mf_dd = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['MF_1','MF_2','MF_3','MF_4','MF_5','MF_6','MF_7','MF_8','MF_9','MF_10','MF_11','MF_12']}
           sum_adult = sum(adult_dd.values())
           sum_juv = sum(juv_dd.values())
           sum_mf = sum(mf_dd.values())
           mf_sum.append(sum_mf)
           
           if (sum_adult == 0) and (sum_juv == 0) and (sum_mf == 0):
               infhost_t += 0
           else:
               infhost_t += 1

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
           meta_popdict[mpop][npop]["A_8"] = random.sample(meta_popdict[mpop][npop]["A_7"],int(round(len(meta_popdict[mpop][npop]["A_7"])*mort_A)))
           meta_popdict[mpop][npop]["A_7"] = random.sample(meta_popdict[mpop][npop]["A_6"],int(round(len(meta_popdict[mpop][npop]["A_6"])*mort_A)))
           meta_popdict[mpop][npop]["A_6"] = random.sample(meta_popdict[mpop][npop]["A_5"],int(round(len(meta_popdict[mpop][npop]["A_5"])*mort_A)))
           meta_popdict[mpop][npop]["A_5"] = random.sample(meta_popdict[mpop][npop]["A_4"],int(round(len(meta_popdict[mpop][npop]["A_4"])*mort_A)))
           meta_popdict[mpop][npop]["A_4"] = random.sample(meta_popdict[mpop][npop]["A_3"],int(round(len(meta_popdict[mpop][npop]["A_3"])*mort_A)))
           meta_popdict[mpop][npop]["A_3"] = random.sample(meta_popdict[mpop][npop]["A_2"],int(round(len(meta_popdict[mpop][npop]["A_2"])*mort_A)))
           meta_popdict[mpop][npop]["A_2"] = random.sample(meta_popdict[mpop][npop]["A_1"],int(round(len(meta_popdict[mpop][npop]["A_1"])*mort_A)))
           # first year adults          
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
               for ad in meta_popdict[mpop][npop]["A_{}".format(i)]:
                   births = np.random.poisson(fecund)                        
                   n = 0
                   while n < births:                         
                       newmf1 = copy.copy(ad)
                       wb_parent2 = meta_popdict[mpop][npop]["A_{}".format(random.randint(1,8))] #random adult                           
                       for p in range(2,len(newmf1)):                             
                           newmf1[p] = newmf1[p] + wb_parent2[p]
                       mf1.append(newmf1)
                       n += 1
           mf = sum(mf1,[])
           
           # recombination in new mf
           num_recomb = []
           for i in range(len(basepairs)):
               num_recomb.append(np.random.binomial(2*len(mf), (basepairs[i]-1) * recombination[i])) 
           if sum(num_recomb) != 0:
               mf = recombination(mf,num_recomb,basepairs)
               
           #makes diploid
           for m in mf:
               for item in m:
                   if len(item) == 4:
                       item = [item[random.randint(0,1)],item[random.randint(2,3)]]
                                    
           #mutation in new mf  
           num_muts = []
           for i in range(0,len(basepairs)):
               if recombination[i] == 0:
                   num_muts.append(np.random.binomial(len(mf), basepairs[i] * mutation[i]))
               else:
                   num_muts.append(np.random.binomial(2*len(mf), basepairs[i] * mutation[i]))               
           if sum(num_muts) != 0:
               mf = mutation(mf,num_muts,basepairs)
           meta_popdict[mpop][npop]["MF_1"] = mf
                
    else: #a year has not passed on months, juveniles and MF move to next age class
       for npop in meta_popdict[mpop].keys(): #inf                
           #count individuals in each class for density dependent calculations
           adult_dd = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['A_1','A_2','A_3','A_4','A_5','A_6','A_7','A_8']}
           juv_dd = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['J_1','J_2','J_3','J_4','J_5','J_6','J_7','J_8','J_9','J_10','J_11','J_12']}
           mf_dd = {key: sum(len(dct[key]) for dct in meta_popdict[mpop][npop].values()) for key in ['MF_1','MF_2','MF_3','MF_4','MF_5','MF_6','MF_7','MF_8','MF_9','MF_10','MF_11','MF_12']}
           sum_adult = sum(adult_dd.values())
           sum_juv = sum(juv_dd.values())
           sum_mf = sum(mf_dd.values())
           mf_sum.append(sum_mf)

           if (sum_adult == 0) and (sum_juv == 0) and (sum_mf == 0):
               infhost_t += 0
           else:
               infhost_t += 1
                      
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
           tempA1 = random.sample(meta_popdict[mpop][npop]["J_12"],int(round(len(meta_popdict[mpop][npop]["J_12"])*mort_J)))
           for subl in tempA1:
               subl[0] += 1
           #first year adults
           meta_popdict[mpop][npop]["A_1"].append(tempA1)              
           # juvenilles                
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

           # birth of MF_1
           mf1 = []                
           for i in range(1,8):
               for ad in meta_popdict[mpop][npop]["A_{}".format(i)]:
                   births = np.random.poisson(fecund)                        
                   n = 0
                   while n < births:                         
                       newmf1 = copy.copy(ad)
                       wb_parent2 = meta_popdict[mpop][npop]["A_{}".format(random.randint(1,8))] #random adult                           
                       for p in range(2,len(newmf1)):                             
                           newmf1[p] = newmf1[p] + wb_parent2[p]
                       mf1.append(newmf1)
                       n += 1
           mf = sum(mf1,[])
           
           # recombination in new mf
           num_recomb = []
           for i in range(len(basepairs)):
               num_recomb.append(np.random.binomial(2*len(mf), (basepairs[i]-1) * recombination[i])) 
           if sum(num_recomb) != 0:
               mf = recombination(mf,num_recomb,basepairs)
               
           #makes diploid
           for m in mf:
               for item in m:
                   if len(item) == 4:
                       item = [item[random.randint(0,1)],item[random.randint(2,3)]]
                                    
           #mutation in new mf  
           num_muts = []
           for i in range(0,len(basepairs)):
               if recombination[i] == 0:
                   num_muts.append(np.random.binomial(len(mf), basepairs[i] * mutation[i]))
               else:
                   num_muts.append(np.random.binomial(2*len(mf), basepairs[i] * mutation[i]))               
           if sum(num_muts) != 0:
               mf = mutation(mf,num_muts,basepairs)
           meta_popdict[mpop][npop]["MF_1"] = mf
            
    return meta_popdict, infhost_t, sum(mf_sum)    

def hostdeath(meta_popdict, mpop):
    '''when a host dies/cures all the MF/Wb disappear from the dict'''
    dpop = random.choice(meta_popdict[mpop].keys()) #draw random host infection
    meta_popdict[mpop][dpop] = {} #blanks the dictionary entry
    return meta_popdict
    
def mutation(mf,num_muts,basepairs):
   '''this is run every time the prob of a mutation is true
   mf:(list) list of mf from adults
   num_muts:(list,int) number of mutations observed  for each locus
   '''      
    
   for i in range(len(basepairs)): #which locus
       muts = 0       #mutation counter
       while muts < num_muts[i]: #keep going until all mutations are assigned
           mut_mf = random.randrange(0,len(mf)) #choose random index in mf list
           new_allele = random.randint(0,basepairs[i])
           if len(mf[mut_mf][i+1]) > 1:
              hap = [random.randint(0,1)]
              if new_allele in mf[mut_mf][i+1][hap]:
                   mf[mut_mf][i+1][hap].remove(new_allele)
              else:
                   mf[mut_mf][i+1][hap].append(new_allele) # add new position
           else:   
              if new_allele in mf[mut_mf][i+1]:
                 mf[mut_mf][i+1].remove(new_allele)
              else:
                 mf[mut_mf][i+1].append(new_allele) # add new position
           muts += 1
   return mf

def recombination(mf,num_recomb,basepairs):
       '''this is run every time the prob of a recombination is true
       num_recomb:(list,int) number of recombination events observed
       mf:(list) list of mf from adults       
       '''      
       for i in range(1,basepairs): #number of muts
          recomb = 0       
          while recomb < num_recomb[i-1]: #keep going until all recombinations are  assigned
               rec_mf = random.randrange(0,len(mf)) #choose random index in mf lis
               new_recomb = random.randint(0,3)
               if new_recomb < 2: #first parent
                      hap1 = mf[rec_mf][i+1][0]
                      hap2 = mf[rec_mf][i+1][1]
                      crossover_pos = random.randint(0,basepairs[i])
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
                      crossover_pos = random.randint(0,basepairs[i])
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
   
def transmission(mpop,transmission_mat,meta_popdict,dispersal,L3trans,hostpopsize):    
    '''transmission function for worms between hosts
    L3trans:(int) number of successful transmission events
    mpop:(int) current metapopulation
    '''
    trans = 0
    mpop = "meta_" + str(mpop+1)
    while trans < L3trans:    
         dpop = random.choice(meta_popdict[mpop].keys()) #random host to donate MF    
         while [len(meta_popdict[mpop][dpop][i]) is 0 for i in meta_popdict[mpop][dpop] if "MF" in i].count(True) is 12: # make ceratin the donating pop actually has MF to donate
             dpop = random.choice(meta_popdict[mpop].keys())    
         if len(transmission_mat[mpop][dpop][1]) >= hostpopsize:
             prob_newinf = 0 #all hosts are infected so no new infections
         else:
             prob_newinf = (float(1)/len(transmission_mat[mpop][dpop][1])) #prob of new infection is determined by number of uninfected  hosts
     
         #new infection
         if random.uniform(0,1) < prob_newinf:
             newpop = "pop_{}".format(len(meta_popdict[mpop].keys()) + 1) #add new infection as population + 1
             meta_popdict[mpop][newpop]["J_1"] = []  #initialize new entry
             rmf = random.randint(1,12) #which age class of MF is the L3 from
             while len(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) is 0:
                 rmf = random.randint(1,12) # keep choosing until MF age pop is not empty
             dmf = random.choice(meta_popdict[mpop][dpop]["MF_{}".format(rmf)])      
             meta_popdict[mpop][newpop]["J_1"].append(dmf)
             transmission_mat = new_infection(transmission_mat,mpop,dpop,newpop,dispersal) #find new transmission position
         #reinfection
         else:     
             rpop = random.choice([i for i, x in enumerate(transmission_mat[mpop][dpop][1]) if x == 1]) + 1 #choose a random value that is 1 as the receiving pop        
             if meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"] is not list:
                 meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"] = []  #in the case that J_1 is empty
             rmf = int(round(random.uniform(.51,12)))
             while len(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) is 0: # in the case that MF is empty in the donating
                 rmf = random.randint(1,12)
             dmf = random.choice(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) #choose random mf      
             meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"].append(dmf) #add dontated mf
         trans += 1
    return meta_popdict, transmission_mat

def new_infection(transmission_mat,mpop,dpop,newpop,dispersal):
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
    return transmission_mat
    
def wb_sims(numberGens,prev,hostpopsize,bites_person,hours2bite,muWormBurden,sizeWormBurden,villages,density_dependence,initial_migration,initial_distance_m,theta,basepairs,mutation,recombination,time2Ancestral,thetaAncestral,thetaRegional,time_join12,time_join23,time_join34,muTrans,sizeTrans,sigma,mortalityHost,mortalityAdult,mortalityJuv,mortalityMF,fecundity):
    '''this will call other functions to intialize, then functions from the life cycle
    numberGens:(int) how long to run the simulation
    burnin:(int) generations after which data is recorded
    prev: (list, float) percent of hosts infected
    hostpopsize:(list, int) number of possible hosts (human population)
    bites_person:(list,int) bites received per person per hour, default = 20
    hours2bite:(list,int) active hours perday for mosquitoes biting, default = 6 ~10PM - 4AM
    **Density dependece for transmission explanation
        Nmf/50 is number in 20ul blood
        culex: L3 = Ks1(1-exp^-(r1*m/Ks1))
        anopheles: L3 = Ks2(1-exp^-(r2(m-T))/Ks2)^2
        where L3 is the number of L3 produced; m is the number of MF per 20ul, r1 is the rate of development with more MF ingested
        Ks1 is the max limiting value of L3 developing.T is the threshold density where above this the MF get a facilitation effect.
        values from gambhir michael 2008: 
             Ks1: 4.406 (+- 0.362)
             r1: 0.019 (+- .058)
             Ks2: 4.395 (+- 0.332)
             r2: 0.055 (+- 0.004)
             T: 0
    '''
    #initialize
    meta_popdict, transmission_mat, dispersal = worm_burden(prev,hostpopsize,muWormBurden,sizeWormBurden,villages,initial_migration,initial_distance_m,theta,basepairs,mutation,recombination,time2Ancestral,thetaAncestral,thetaRegional,time_join12,time_join23,time_join34,muTrans,sizeTrans,sigma)  
    #set counters
    time_month = 0 #begin time
    sim_time = numberGens #gens to pass
    infhost = [i*j for i,j in zip(prev,hostpopsize)] #staring prevalance
    prev_t = prev
    mf_mpopavg = [0] * villages #0s for MF
    while time_month <= sim_time:
        for mpop in range(villages):
             #maturation
             meta_popdict, infhost[mpop], mf_mpopavg[mpop] = maturation(mpop,meta_popdict,time_month,infhost[mpop],density_dependence,mortalityHost,mortalityAdult,mortalityJuv,mortalityMF,fecundity,basepairs,mutation,recombination) #update matrices   
             #transmission
             totalbitesvillage = bites_person[mpop] * hours2bite[mpop] * 30 * hostpopsize[mpop] #20 bites per person per hour * 6 hours (atnight) * 30 days * hostpopsize
             infbites=np.random.binomial(totalbitesvillage,(prev_t[mpop]*0.37)) #prob of bite on infected person picking up MF. prob of 0.37 from Gambhir paper
             if density_dependence: #values for anopheles
                  mf_blood = (((mf_mpopavg[mpop])/infhost[mpop])/235)/50 #(sum[allMFstages]/235)/50 number of MF in 20ul of blood; 235ml is 5% of total body blood and there are 50*20ul in 1 ML
                  L3trans = round(infbites*(4.395*(1-math.exp(-(0.055*(mf_blood))/4.395))**2) * (0.414*0.32)) #proportion of L3 that leave mosquito is 0.414 and prop that enter host after leaving mosquito is 0.32
             else:
                  L3trans = np.random.binomial(infbites,(0.414*0.32)) #prob of infected bite that picks up MF which then mature to L3 and spread to another host
             meta_popdict, transmission_mat = transmission(mpop,transmission_mat,meta_popdict,dispersal,L3trans,hostpopsize[mpop])  
        prev_t[mpop] = infhost[mpop]/hostpopsize[mpop]
        time_month += 1
        
#def main():
#    wb_sims(args.*)
#if __name__ == '__main__':
#    main()

'''unused gillespie algorithm for continuous transmission:
   #random_t = np.random.uniform(0,1,len(rates)) #vector of random variables
   #wait_t = -(np.log(random_t))/rates #wait time to next event
   #event_t = np.argmin(wait_t) #which event has the shortest wait time
   #time =+ wait_t[event_t] #update time
   #rate of bite is: 20 bites per person per month * hostpopsize / 720 = bites/hour; rate_b = (20.00*hostpopsize[0])/720   
   #rate of bites that have MF: rate of bite * prev * prob pickup MF; rate_ib = rate_b * len(meta_popdict.keys()) * 0.37 
   #rate of transmission: rate_ib * number_L3 * reinjected to new host; number_L3 = 4.395(1-exp^-(.055(m))/4.395)^2; leave with bite = .414; reinjected = .32; m=(sum[donor["MF"]]/235)/50    
   #(sum[allMFstages]/235)/50 number of MF in 20ul of blood 235ml is 5% of total body blood and there are 50*20ul in 1 ML
    '''


