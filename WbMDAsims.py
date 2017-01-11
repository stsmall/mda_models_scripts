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
    hostdeath(meta_popdict, *mpop)
    recombination(*mf, *num_recomb, basepairs)
    mutation(*mf, *num_mut, basepairs)
  transmission(transmission_mat, *mpop, meta_popdict, dispersal)
    new_infection(transmission_mat, mpop, dpop, newpop, dispersal)
##
@author:stsmall
"""
import math
import random
import argparse
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import time
from Ipython import embed
import pandas as pd

from funcs import ms_outcall
from recombination import recombination_fx
from maturation import maturation


#meta_popdict, transmission_mat, dispersal, hap_pop = worm_burden([10,10], [5,5], [50,5], [[1,1],[2,2]], [1000,2000], [7.6E-8, 2.9E-9],  [0, 2.9E-9], 1800, 344, 23, 240, 240, 240, 2, 0.0001, [1000], 100, 1, 2)
#meta_popdict, X, Y = maturation(1, meta_popdict, 10, True, .90, .90, .90, .90, 20, [1000, 2000], [7.6E-8, 2.9E-9],  [0.0, 2.9E-9])

def get_args():
    '''asdfasdf'''
    parser = argparse.ArgumentParser()
    ## initialize
    #migration_matrix
    parser.add_argument('-im', '--initial_migration', type=float, default=.0001, 
            help=("migration rate between villages/metapopulations in the" 
            "model.This is strictly for initial conditions and can be" 
            "changed for the forward-in-time portion"))
    parser.add_argument('-idm', '--initial_distance_m', type=list, 
            help="initial_distance_m is list [1000] such that distance_m[0] is between 1 & 2")
    parser.add_argument('-v', '--villages', type=int, default=1, 
            help="sets the intial number of villages.")
    parser.add_argument('-t', '--theta', type=list, required=True, 
            help=("observed theta value of worm populations for" 
            "each locus in format [[meta1_locus1,meta2_locus1],"
            "[meta1_locus2, meta2_locus2]]"))
    parser.add_argument('-bp', '--basepairs', type=list, default=13000, 
            help="length in basepairs of each locus")
    parser.add_argument('-u', '--mutation', type=list, default=7.6E-8, 
            help="expected as list, mutation rate per bp per generation for each locus")
    #ms_outcall
    parser.add_argument('-t12', '--time_join12', type=int, default=240, 
            help="generations until time of joining for pop1 and pop2")
    parser.add_argument('-t23', '--time_join23', type=int, default=240, 
            help="generations until time of joining for pop2 and pop3")
    parser.add_argument('-t34', '--time_join34', type=int, default=240, 
            help="generations until time of joining for pop3 and pop4")
    parser.add_argument('-t2a', '--time2Ancestral', type=int, default=1800, 
            help="generations until ancestral population for PNG is 1800 generations for Africa/Haiti 500 generations")
    parser.add_argument('-at', '--thetaAncestral', type=int, default=344, 
            help="ancestral theta before")
    parser.add_argument('-ar', '--thetaRegional', type=int, default=23, 
            help="regional theta")
    parser.add_argument('-r', '--recombination', type=list, default=0, 
            help="recombination for each locus, if 0 assumes haploid is 1")
    #trans_init
    parser.add_argument('-mt', '--muTrans', type=int, default=100, help="mu for neg bino in transmission, distances between hosts")
    parser.add_argument('-st', '--sizeTrans', type=int, default=1, help="size for neg bino in transmission, distance between hosts")
    parser.add_argument('-dp', '--sigma', type=int, default=2, help="sigma, for dispersal")
    #worm_burden
    parser.add_argument('-prev', '--prev', type=list, required=True, help="prevalance of Wb in host populations; should be a list")
    parser.add_argument('-host', '--hostpopsize', type=list, required=True, help="size of host population")
    parser.add_argument('-mw', '--muWormBurden', type=list, help="mu for negative binomial in worm burden per host, add value or defined by uniform")
    parser.add_argument('-sw', '--sizeWormBurden', type=list, help="size for negative binomial in worm burden per host")
    ## wb_sims
    parser.add_argument('-ng', '--numberGens', type=int, required=True, help="total number of generation to run the simulation")
    parser.add_argument('-bpp', '--bites_person', type=int, default=20, help="number of bites recieved per person per hour")
    parser.add_argument('-h2b', '--hours2bite', type=int, default=6, help="number of hours exposed to mosquitoes")
    #wblifecycles
    parser.add_argument('-ma', '--mortalityAdult', type=float, default=0.8863848717161292, help="adult worms dies at rate 0.01 per month or survive at .99")
    parser.add_argument('-mj', '--mortalityJuv', type=float, default=0.865880173963825, help="prob that juv makes it to adult is 0.2 or 0.8177651 per month")
    parser.add_argument('-mf', '--mortalityMF', type=float, default=0.90, help="MF die at rate 0.10 per month or survive at 0.90")
    parser.add_argument('-mh', '--mortalityHost', type=float, default=0.014286, help="Host death per year")
    parser.add_argument('-f', '--fecundity', type=int, default=20, help="mean number of MF born to each Adult per month")
    parser.add_argument('-Dd', '--density_dependence', action="store_true", help="use density dependence")
    #not used
    parser.add_argument('-gtime', '--generations', type=int, default=0.125, help="generation time in years")
    parser.add_argument('-hostmig', '--host_migration_rates', help="list of host migration rates between villages per month")
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

    

def transmission_init(infhost, muTrans, sigma, sizeTrans):
    '''Creates a transmission matrix for locations of infected hosts

    Parameters
    ----------
    infhost: 
    muTrans:(float) 
        parameter of neg binomial
    sizeTrans:(float) 
        parameter of neg binomial
    sigma:(float) 
        parameter for calculating dispersal distance

    Returns
    -------
    transmission_matrix: dataframe
        polar coordinate positions for pops
    dispersal: float
    '''
    print("starting trans_init")
    #set dispersal parameters
    dispersal = 4*math.pi*sigma**2
    pos = np.random.negative_binomial(sizeTrans,
            sizeTrans/float((sizeTrans+muTrans), (infhost, 2)))
    # :TODO fast polar conversion
    pospolar = to_polar(pos)

    transmission_matrix = pd.DataFrame({
                                        'hostpop' : range(infhost),
                                        'angle': pospolar[:, 0], 
                                        'radius': pospolar[:, 1],},
                                       )
    #{'meta_1': {'pop_1': ['x,y',[0,1,1,1,1,1]],'pop_2': ['x,y',[0,1,1,1,0,0]]}}
    return transmission_matrix, dispersal



def worm_burden(infhost, muWormBurden, sizeWormBurden, theta, basepairs, 
        mutation, recombination, time2Ancestral, thetaAncestral, 
        thetaRegional, time_join12, time_join23, time_join34, villages, 
        initial_migration, initial_distance_m, muTrans, sizeTrans, sigma):
    """Worm_burden calls all the previous functions to create meta_popdict.

    Parameters
    ----------
    meta_popdict has the age/stage structure as well as the haplotypes for each locus.
    transmission_mat gives the corresponding locations of the infected hosts
    infhost : (list, int) number of possible hosts (human population)
    muWormBurden: (list,int) 
        avg_burden, average number of adult female worms in infections;
    sizeWormBurden:(list,int) 
        dispersion, size parameter for negative binomial distribution

    Returns
    -------
    """
    print("starting worm_burden")
    t0 = time.clock()
    # worm burden (number of adult worms) per host
    pop_init = []
    for mu, size, numworms in zip(muWormBurden, sizeWormBurden, infhost):
        wb_burden = np.random.negative_binomial(mu, mu/float(mu+size), numworms) # number of successes, prob of success (size/size+mu),number of values to return
        pop_init.append(np.array(wb_burden).tolist())

    # call to function msoutcall to create hap_pop list
    worm_popsize = []
    for meta in pop_init:
        worm_popsize.append(sum(meta))
    hap_pop = ms_outcall(worm_popsize, villages, initial_migration, initial_distance_m, theta, basepairs, mutation, recombination, time2Ancestral, thetaAncestral, thetaRegional, time_join12, time_join23, time_join34)

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
            meta_popdict["meta_" + str(meta_n)]["pop_" + str(pop)]["A_1"] = []  #only A_1 worms
            while j < hostpop:
                haplotypes = [np.random.uniform(), [hap_pop["locus0"][k]]] #initializes info for each worm [rand_num,[],[,],[,],[,],[,]]
                for i in range(1, locus):
                    haplotypes.append([hap_pop["locus"+str(i)][kd], hap_pop["locus"+str(i)][kd+1]])
                meta_popdict["meta_" + str(meta_n)]["pop_" + str(pop)]["A_1"].append(haplotypes)
                j += 1 #counts the A1 in pop_init
                k += 1 #counts the haps in hap_pop[]
                kd += 2 #counts the diploid haps in hap_pop[]
            pop += 1 #advances the pop counter
        meta_n += 1 #advances the meta counter
        pop = 1 #resets pop counter for new meta populations aka village
    transmission_mat, dispersal = transmission_init(infhost, muTrans, sigma, sizeTrans)
    print(time.clock()-t0)
    return meta_popdict, transmission_mat, dispersal, hap_pop
    #meta_1:pop1:[rand,[locus0],[locus1],[locusN]]
    #{'meta_1': {'pop_1': {'A_1': [rand, [[L1_hap]], [[L2_hap1], [L2_hap2]], [[L3_hap1], [L3_hap2]]], [rand, [L1_hap], [L2_hap1, L2_hap2], [L3_hap1, L3_hap2]]

############################################
##Begin simulation functions
############################################



def mutation_fx(mf, num_muts, basepairs):
    '''this is run every time the prob of a mutation is true
   mf:(list) list of mf from adults
   num_muts:(list,int) number of mutations observed  for each locus
   '''
    for i, bp in enumerate(basepairs): #which locus
        muts = 0       #mutation counter
        while muts < num_muts[i]: #keep going until all mutations are assigned
            mut_mf = random.randrange(0, len(mf)) #choose random index in mf list
            new_allele = random.randint(0, bp)
            if len(mf[mut_mf][i+1]) > 1:
                hap = [random.randint(0, 1)]
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


def transmission_fx(mpop, transmission_mat, meta_popdict, dispersal, L3trans, hostpopsize):
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
        if random.uniform(0, 1) < prob_newinf:
            newpop = "pop_{}".format(len(meta_popdict[mpop].keys()) + 1) #add new infection as population + 1
            meta_popdict[mpop][newpop]["J_1"] = []  #initialize new entry
            rmf = random.randint(1, 12) #which age class of MF is the L3 from
            while len(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) is 0:
                rmf = random.randint(1, 12) # keep choosing until MF age pop is not empty
            dmf = random.choice(meta_popdict[mpop][dpop]["MF_{}".format(rmf)])
            meta_popdict[mpop][newpop]["J_1"].append(dmf)
            transmission_mat = new_infection(transmission_mat, mpop, dpop, newpop, dispersal) #find new transmission position
        #reinfection
        else:
            rpop = random.choice([i for i, x in enumerate(transmission_mat[mpop][dpop][1]) if x == 1]) + 1 #choose a random value that is 1 as the receiving pop
            if meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"] is not list:
                meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"] = []  #in the case that J_1 is empty
            rmf = random.randint(1, 12)
            while len(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) is 0: # in the case that MF is empty in the donating
                rmf = random.randint(1, 12)
            dmf = random.choice(meta_popdict[mpop][dpop]["MF_{}".format(rmf)]) #choose random mf
            meta_popdict[mpop]["pop_{}".format(rpop)]["J_1"].append(dmf) #add dontated mf
        trans += 1
    return meta_popdict, transmission_mat

def new_infection(transmission_mat, village, dpop, newpop, dispersal):
    '''
    Parameters
    ----------
    transmission_mat: np.array
        pairwise eucileadean distance between people 
    village: int
        metapopulation index
    dpop: int
        donation population, where genotypes are coming
    newpop: strings
        receiving population
    this is run every time a new individual is infected to rebuild the transmission prob matrix, updates'''
    x1n = np.random.uniform(-1*math.sqrt(dispersal), math.sqrt(dispersal))
    x2n = np.random.uniform(-1*math.sqrt(dispersal), math.sqrt(dispersal))
    # :TODO transmission matrix to dataframe

    transmission_matrix = pd.DataFrame({'village' : [],
                                        'hostpop' : [],
                                        'angle': [], 
                                        'radius': [],},
                                        index = ['1_p1'])

    distances = transmission_matrix.ix[x, y].distance()




    transmission_mat[village][newpop] =\
    ["{},{}".format(int(transmission_mat[village][dpop][0].split(",")[0]) + x1n, int(transmission_mat[mpop][dpop][0].split(",")[1]) + x2n), []]
    for tpop in range(1, len(transmission_mat[mpop].keys())):
        dist = math.sqrt((float(transmission_mat[mpop]["pop_{}".format(tpop)][0].split(",")[0]) - x1n)**2 + (float(transmission_mat[mpop]["pop_{}".format(tpop)][0].split(",")[1]) - x2n)**2)
        if dist < dispersal:
            dist = 1
        else:
            dist = 0
        transmission_mat[mpop][newpop][1].append(dist)
        transmission_mat[mpop]["pop_{}".format(tpop)][1].append(dist)
    transmission_mat[mpop][newpop][1].append(1)
    return transmission_mat

def wb_sims(numberGens, prev, hostpopsize, villages, bites_person, hours2bite, muWormBurden, sizeWormBurden, theta, basepairs, mutation, recombination, time2Ancestral, thetaAncestral, thetaRegional, time_join12, time_join23, time_join34, initial_migration, initial_distance_m, muTrans, sizeTrans, sigma, density_dependence, mortalityHost, mortalityAdult, mortalityJuv, mortalityMF, fecundity):
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
    #set counters
    time_month = 0 #begin time
    sim_time = numberGens #gens to pass
    infhost = [round(i*j) for i, j in zip(prev, hostpopsize)] #staring prevalance
    prev_t = prev
    mf_mpopavg = [0] * villages #0s for MF
    #initialize
    meta_popdict, transmission_mat, dispersal = worm_burden(infhost, muWormBurden, sizeWormBurden, villages, initial_migration, initial_distance_m, theta, basepairs, mutation, recombination, time2Ancestral, thetaAncestral, thetaRegional, time_join12, time_join23, time_join34, muTrans, sizeTrans, sigma)
    while time_month <= sim_time:
        for mpop in range(villages):
             #maturation
            meta_popdict, inf, mfsum = maturation(mpop, meta_popdict, time_month, infhost[mpop], density_dependence, mortalityHost, mortalityAdult, mortalityJuv, mortalityMF, fecundity, basepairs, mutation, recombination) #update matrices
            infhost[mpop] = inf
            mf_mpopavg[mpop] = mfsum

            #transmission
            totalbitesvillage = bites_person[mpop] * hours2bite[mpop] * 30 * hostpopsize[mpop] #20 bites per person per hour * 6 hours (atnight) * 30 days * hostpopsize
            infbites = np.random.binomial(totalbitesvillage, (prev_t[mpop]*0.37)) #prob of bite on infected person picking up MF. prob of 0.37 from Gambhir paper
            if density_dependence: #values for anopheles
                mf_blood = (((mf_mpopavg[mpop])/infhost[mpop])/235)/50 #(sum[allMFstages]/235)/50 number of MF in 20ul of blood; 235ml is 5% of total body blood and there are 50*20ul in 1 ML
                L3trans = round(infbites*(4.395*(1-math.exp(-(0.055*(mf_blood))/4.395))**2) * (0.414*0.32)) #proportion of L3 that leave mosquito is 0.414 and prop that enter host after leaving mosquito is 0.32
            else:
                L3trans = np.random.binomial(infbites, (0.414*0.32)) #prob of infected bite that picks up MF which then mature to L3 and spread to another host

            meta_popdict, transmission_mat = transmission_fx(mpop, transmission_mat, meta_popdict, dispersal, L3trans, hostpopsize[mpop])
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
embed()
