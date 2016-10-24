# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:56:02 2016
A simulation to examine population genetic statistics in a epidemiological setting
incorporating historical demographic history in the parasite to account for observed patterns of genetic diversity.
There is a ridiculous number of options to set so see --help.
NOTE: If you want to use the gillepsie algorithm invoked by --continuous then it can only be 
a sinlge village since I couldnt figure out how to sync the times with waiting times. 
For multiple villages time is in discrete units
@author: stsmall
"""
#install anaconda
import numpy as np
from math import pi, exp, sqrt
import subprocess, re, random
from sklearn.metrics.pairwise import euclidean_distances
import argparse
import copy

#-----------------------------------------------------------------------------
# arguments
def get_args():
    parser = argparse.ArgumentParser()
    #class wbsims
    parser.add_argument('-ng','--numberGens', type=int, required=True, help="total number of generation to run the simulation")
    #class wbinit
    parser.add_argument('-p', '--prev', type=list, required=True, help="prevalance of Wb in host populations; should be a list")
    parser.add_argument('-h', '--hostpopsize',type=list, required=True,help="size of host population")
    parser.add_argument('-d', '--distance',type=list, required=True, help="distance between host populations")
    parser.add_argument('-t', '--theta',type=list,required=True, help="theta value of worm populations")
    parser.add_argument('-a','--thetaAncestral', type=int,default=10, help="ancestral theta before split of worm populations")
    parser.add_argument('-mw', '--muWormBurden',type=int, help="mu for negative binomial in worm burden per host, add value or defined by uniform")
    parser.add_argument('-sw', '--sizeWormBurden',type=int, help="size for negative binomial in worm burden per host")
    parser.add_argument('-mt', '--muTrans',type=int,default=100, help="mu for neg bino in transmission, distances between hosts")
    parser.add_argument('-st', '--sizeTrans',type=int,default=1, help="size for neg bino in transmission, distance between hosts")
    parser.add_argument('-dp', '--sigma',type=int,default=2, help="sigma, for dispersal")
    parser.add_argument('-t12', '--time12',type=int,default=0.05, help="time of joining for pop1 and pop2")
    parser.add_argument('-t23', '--time23',type=int,default=0.05,help="time of joining for pop2 and pop3")
    parser.add_argument('-t45', '--time34',type=int,default=0.05,help="time of joining for pop3 and pop4")
    parser.add_argument('-ta','--timeA',type=int, default=0.05)
    parser.add_argument('-m', '--migration',type=int,default=.0001,help="migration rate between hostpops")
    parser.add_argument('-N', '--Ne',type=int,default=10000,help="effective population size")
    #class wblifecycles
    parser.add_argument('-ma','--mortalityAdults',type=int,default=0.8863848717161292,help="adult worms dies at rate 0.01 per month or survive at .99")
    parser.add_argument('-mj','--mortalityJuv',type=int,default=0.865880173963825,help="prob that juv makes it to adult is 0.2 or 0.8177651 per month")
    parser.add_argument('-mf','--mortalityMF',type=int,default=0.90,help="MF die at rate 0.10 per month or survive at 0.90")
    parser.add_argument('-f','--fecundity',type=int,default=20,help="")
    parser.add_argument('-u','--mutation',type=int,default=7.6E-8,help="")
    parser.add_argument('-g','--generations',type=int,default=0.125,help="generation time in years")
    parser.add_argument('-bp','--basepairs',type=int,default=13000,help="")
    parser.add_arguement('-cont','--continuous',action='store_true',help='runs transmission as continuous using gillespie algorithm and state dependent rates. This only works for 1 village and invoking will default everything to 1 village')
    parser.add_arguement('-hostmig','--host_migrationrates', help="list of host migration rates between villages per month")
    args = parser.parse_args()
    return args

#------------------------------------------------------------------------------
# class for dictionaries

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
            
#-----------------------------------------------------------------------------
# functions for model intialization           
            
class wbinit:
    def __init__(self,args):
        self.prev = args.prev
        self.hostpopsize = args.hostpopsize
        self.distance = args.distance
        self.muWormBurden = args.muWormBurden
        self.sizeWormBurden = args.sizeWormBurden
        self.theta = args.theta
        self.thetaAncestral = args.thetaAncestral
        self.sigma = args.sigma
        self.muTrans = args.muTrans
        self.sizeTrans = args.sizeTrans
        self.timeA = args.timeA
        self.time12 = args.time12
        self.time23 = args.time23
        self.time34 = args.time34
        self.migration = args.migration
        self.Ne = args.Ne
        self.continuous = args.continuous
        self.hostmig = args.host_migrationrates

    def migration_matrix(self,villages):
        '''calculated a migration matrix for use with ms (hudson 2000) from euclidian
        distances. The value of 2Nm is weighted by the distance from the next population
        as an exponential random variable. The highest this can be is 2Nm/1'''  
        #distance_m is list [1000,1500] such that distance_m[0] is between 1 & 2
        #distance_m[1] is between 1 & 3 etc ...
        #villages is int 3; such that len(distance_m) == ((villages)*(villages-1)/2)
        m = self.migration
        Ne = self.Ne
        distance_m = self.distance
        if villages > 4: 
            raise ValueError("only handles 4 villages ATM")
        elif villages < 4:
            if len(distance_m) != ((villages)*(villages-1)/2): raise ValueError("there are not adequate pairwise comparisons")
            mig = []
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
    
    def ms_outcall(self,worm_popsize):
        '''external call to ms (Hudson 2000) using the migration matrix calculated
    above from the euclidian distances. This function will then read in the stdout from
    ms. The continious location of Seg given by ms is transformed to discrete coordinates. 
    It then builds a list of lists from the discrete coordinates which is then pointer
    referenced in the dictionary recording the metapopulation information
        metapop_init = [100,100] #2 metapops
        theta = [metapop1,metapop2]    
        ''' 
        #to add multiple loci just change howmany to >2 this will produce a few poistionis lines which then make multiple hap_pop dicts and have meta_popdict point to each haplotype
        #in mutation you will need to change the rates and rerun mutation for each locus, namely each hap_pop dictionary
        theta = self.theta
        theta_anc = self.thetaAncestral
        ta = self.timeA
        t12 = self.time12
        t23 = self.time23
        t34 = self.time34
        if len(worm_popsize) == 1: #1 village
            mscmd = "ms {} 1 -t {} -eN {} {} > temp_file".format(worm_popsize[0],theta[0],ta,float(theta_anc/theta[0]))        
        else: #set up for 2 villages
            num_subpops = len(worm_popsize) #-I num_pops
            total_inds = sum(worm_popsize) # how many
            sub_pop = " ".join(map(str,worm_popsize))#-I X i j ...
            #t = 0.05 #2000gens/4*N0  #merge time   
            #-ej the villages split at 2000 gens in the past which is 2000 years
            if len(worm_popsize) == 2:        
                mscmd = "ms {} 1 -t {} -I {} {} -n 1 {} -n 2 {} -ma {} -ej {} 1 2 -eN {} {} > temp_file".format(total_inds,theta[0],num_subpops,sub_pop,1,float(theta[1])/theta[0],self.migration_matrix(len(worm_popsize)),t12,ta,theta_anc)       
            #after the split the pops are half the size of ancestral and share migrants at rate -ma
            #pop continues to go back until ancestral time    
            elif len(worm_popsize)==3:
                mscmd = "ms {} 1 -t {} -I {} {} -n 1 {} -n 2 {} -n 3 {} -ma {} -ej {} 1 2 -ej {} 2 3 -eN {} 1 > temp_file".format(total_inds,theta[0],num_subpops,sub_pop,1,float(theta[1])/theta[0],float(theta[2])/theta[0],self.migration_matrix(len(worm_popsize)),t12,t23,ta,theta_anc)       
            elif len(worm_popsize)==4:
                mscmd = "ms {} 1 -t {} -I {} {} -n 1 {} -n 2 {} -n 3 {} -n4 {} -ma {} -ej {} 1 2 -ej {} 2 3 -ej {} 3 4 -eN {} 1> temp_file".format(total_inds,theta[0],num_subpops,sub_pop,1,float(theta[1])/theta[0],float(theta[2])/theta[0],float(theta[3])/theta[0],self.migration_matrix(len(worm_popsize)),t12,t23,t34,ta,theta_anc)       
        print mscmd
        proc = subprocess.Popen(mscmd, shell=True)
        proc.wait()
    
        #parses ms output    
        hap_pop = []  
        with open("temp_file",'r') as ms:
            for line in ms:
                if line.startswith("positions"):             
                    for line in ms:            
                        hap = line.rstrip("\n")                    
                        hap_pop.append([m.start() for m in re.finditer("1", hap)])    
        return hap_pop
        
    def trans_init(self):
        '''initializes locations for above infections within
        prevalence = [.8,.2]
        hostpop_size = [100,100]    
        '''
        #from sklearn.metrics.pairwise import euclidean_distances
        sigma = self.sigma
        dispersal = 4*pi*sigma**2
        if self.sizeTrans:
            size = self.sizeTrans
        else:
            size = 1
        if self.muTrans:
            mu = self.muTrans
        else:
            mu = 100
        metan = 1
        transmission_mat = AutoVivification()     
        for i,j in zip(self.prev,self.hostpopsize):    
            x1 = np.random.negative_binomial(size,size/float((size+mu)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
            x2 = np.random.negative_binomial(size,size/float((size+mu)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
            X = np.vstack((x1,x2)).T    
            dist = euclidean_distances(X,X) 
            for pop in range(0,len(x1)):
                transmission_mat["meta_{}".format(metan)]["pop_{}".format(pop+1)] = ["{},{}".format(x1[pop],x2[pop]), [1 if item < dispersal else 0 for item in dist[pop,:]]]
            metan += 1
        return transmission_mat, dispersal    
      
    def worm_burden(self):
        '''
        num_pops = [15,25] #number of infections in meta1 and meta2; prevalance * population
        mu = [5,1] #avg_burden, average number of adult female worms in infections;
        size = [1,50] #size = dispersion, size parameter for negative binomial distribution
        hap_pop is from ms_outcall    
        '''  
        # number of adult worms per infection    
        num_inf = []
        for i,j in zip(self.prev,self.hostpopsize):
            num_inf.append(round(i*j))
        pops = len(self.hostpopsize)
        if self.muWormBurden:
            mu = self.muWormBurden            
        else:
            mu = np.random.uniform(1,5,pops) #return num equal to
        if self.sizeWormBurden:
            size = self.sizeWormBurden
        else:
            size = np.random.uniform(1,50,pops) #return num equal to
        #worm burden    
        pop_init=[]
        for i,j,k in zip(mu,size,num_inf):
            wb_burden = np.random.negative_binomial(i,i/float(i+j),k) # number of successes, prob of success (size/size+mu),number of values to return
            pop_init.append(np.array(wb_burden).tolist())
        
        #make populations with haplotypes
        worm_popsize = []
        for meta in pop_init:    
            worm_popsize.append(sum(meta))
        hap_pop = self.ms_outcall(worm_popsize)  
            
       #meta_popdict[meta1][pop][age][list]
        meta_popdict = AutoVivification()            
        pop = 1 #initial value
        meta = 1 #initial value
        k = 0 #initial start of happop
        for metapop in pop_init:
            for wb_a1 in metapop:
                j = 0 #initializes         
                meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"]=[]
                while j < wb_a1:
                    meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"].append([np.random.uniform(), k])
                    j += 1 #counts the A1 in population
                    k += 1 #counts the haps in hap_pop[]
                pop += 1 #advances the pop counter
            meta += 1 #advances the meta counter
            pop = 1 #resets pop counter for new meta
        transmission_mat,dispersal=self.trans_init(self.prev,self.hostpopsize)
        return meta_popdict, hap_pop, transmission_mat #dictionary, max(max(hap_pop))

#----------------------------------------------------------------------------
# life cycle

class wblifecycle:
    def __init__(self,args):
        self.mortalityAdults = self.mortalityAdults
        self.mortalityJuv = self.mortalityJuv
        self.mortalityMF = self.mortalityMF
        self.fecundity = self.fecundity
        self.mutation = self.mutation
        self.generations = self.generations
        self.basepairs = self.basepairs
    
    def maturation(self, meta_popdict, hap_pop, month):
        '''step all individuals forward, after deaths. Kill them then move, then births'''    

        # density dependent mortality    

    #    b = 1.0/K #K is carrying capacity
    #    S = 0.99
    #    a = (S*b)/exp(-1) #where S is maximum survival
    #    mort_A = a * sum(len(A_1 ... A_8)) * exp(-b * sum(len(A_1 ... A_8))) #Ricker fx 

    # fixed mortality    

        mort_A = self.mortalityAdults #year
        mort_J = self.mortalityJuv #month
        mort_M = self.mortalityMF #month
        
    # density dep fecundity   

    #    b = 1.0/K #K is carrying capacity
    #    S = 0.99
    #    a = (S*b)/exp(-1) #where S is maximum survival
    #    fecund = 20 * a * sum(len(A_1 ... A_8)) * exp(-b * sum(len(A_1 ... A_8))) #Ricker fx 

    # fixed fecundity     

        fecund = self.fecundity     
    
    # mutation

        bp = self.basepairs
        pmut = self.mutation
        frac_gen = self.generations #gen is 8 months so 1 month is .125 of a generation
        
        rzero_freq = []
        hap_freq = []
       #since this month to month 
        if month%12 is 0: #this denotes 1 year has passed so adults mature to next age class 
            for mpop in meta_popdict.keys(): #villages
                for npop in meta_popdict[mpop].keys(): #inf
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
                    #count A_1                
                    rzero_freq.append([i[0]for i in meta_popdict[mpop][npop]["A_1"]])
                    #reset A_1 since they are new from Juv12
                    for subl in meta_popdict[mpop][npop]["A_1"]:
                        subl[0] = random.random()
                    #MF_1
                    mf1 = []                
                    for i in range(1,8):
                        for j in meta_popdict[mpop][npop]["A_{}".format(i)]:
                            mf1.append([j]*np.random.poisson(fecund)) 
                    mf = sum(mf1,[])
                    #set up mutation        
                    num_muts = np.random.binomial(len(mf), bp * pmut * frac_gen)
                    if num_muts != 0:
                        mf,hap_pop = self.mutation(mf,hap_pop,num_muts)
                    meta_popdict[mpop][npop]["MF_1"] = mf
                    hap_freq.append([i[1]for i in mf])
                    
        else: #a year has not passed on months, juveniles and MF move to next age class
            for mpop in meta_popdict.keys():        
                for npop in meta_popdict[mpop].keys(): #inf
                    # first year adults, later add from different age classes with maturity an evolvable trait            
                    meta_popdict[mpop][npop]["A_1"].append(random.sample(meta_popdict[mpop][npop]["J_12"],int(round(len(meta_popdict[mpop][npop]["J_12"])*mort_J))))
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
                    #count A_1                
                    rzero_freq.append([i[0]for i in meta_popdict[mpop][npop]["A_1"]])
                    #reset A_1
                    for subl in meta_popdict[mpop][npop]["A_1"]:
                        subl[0] = random.random()
                    # MF_1
                    mf1 = []                
                    for i in range(1,8):
                        for j in meta_popdict[mpop][npop]["A_{}".format(i)]:
                            mf1.append([j]*np.random.poisson(fecund)) 
                    mf = sum(mf1,[])
                    #set up mutation        
                    num_muts = np.random.binomial(len(mf), bp * pmut * frac_gen)
                    if num_muts != 0:
                        mf,hap_pop = self.mutation(mf,hap_pop,num_muts)
                    meta_popdict[mpop][npop]["MF_1"] = mf
                    hap_freq.append([i[1]for i in mf])
    
        #this will calculate for the enitre meta and each inf
        return meta_popdict, hap_pop, {x:rzero_freq.count(x) for x in rzero_freq},{y:hap_freq.count(y) for y in hap_freq}  
        
    def mutation(self,mf,hap_pop,num_muts):
       '''this is run every time the prob of a mutation is true, updates seq_base'''   
       #calculate the number of mutations expected
       mut_mf = [random.randrange(len(mf)) for i in range(num_muts)] #choose random index in mf list
           #add new sequence to hap_pop
       for m in mut_mf:                 
           new_hap = copy.copy(hap_pop[mf[m][1]])
           new_allele = (max([max(a) for a in hap_pop])) + 1
           new_hap.append(new_allele)
           hap_pop.append(new_hap)
           mf[m][1] = len(hap_pop)-1
       return mf, hap_pop
    
    def transmission(self,transmission_mat,meta_p,meta_popdict,dispersal):    
        '''this has a continious prob using a gillepsie algorithm to get the waiting time between transmission events'''    
        #this executes if tranmission is true     
        #pick a random donating pop:
        dpop = random.choice(meta_popdict[meta_p].keys())    
        while [len(meta_popdict[meta_p][dpop][i]) is 0 for i in meta_popdict[meta_p][dpop] if "MF" in i].count(True) is 12:
            dpop = random.choice(meta_popdict[meta_p].keys())
        r = random.uniform(0,1)    
        #new infection
        if r > 1 - (float(1)/len(transmission_mat[meta_p][dpop][1])):
            newpop = "pop_{}".format(len(meta_popdict[meta_p].keys())+1)
            meta_popdict[meta_p][newpop]["J_1"] = []
            rmf = int(round(random.uniform(.51,12)))
            while len(meta_popdict[meta_p][dpop]["MF_{}".format(rmf)]) is 0:
                rmf = int(round(random.uniform(.51,12)))
            dmf = random.choice(meta_popdict[meta_p][dpop]["MF_{}".format(rmf)])      
            meta_popdict[meta_p][newpop]["J_1"].append(dmf)
            print "new"
            transmission_mat = self.new_infection(transmission_mat, meta_p,dpop,newpop,dispersal) #find new transmission position
        #reinfection
        else:     
            rpop = random.choice([i for i, x in enumerate(transmission_mat[meta_p][dpop][1]) if x == 1]) + 1 #choose a random value that is 1 as the receiving pop        
            if meta_popdict[meta_p]["pop_{}".format(rpop)]["J_1"] is not list:
                meta_popdict[meta_p]["pop_{}".format(rpop)]["J_1"] = []
            rmf = int(round(random.uniform(.51,12)))
            while len(meta_popdict[meta_p][dpop]["MF_{}".format(rmf)]) is 0:
                rmf = int(round(random.uniform(.51,12)))
            dmf = random.choice(meta_popdict[meta_p][dpop]["MF_{}".format(rmf)])      
            meta_popdict[meta_p]["pop_{}".format(rpop)]["J_1"].append(dmf)
        return meta_popdict, transmission_mat
    
    def new_infection(self,transmission_mat, meta_p,dpop,newpop,dispersal):
        '''this is run every time a new individual is infected to rebuild the transmission prob matrix, updates'''        
        x1n = np.random.uniform(-1*sqrt(dispersal), sqrt(dispersal)) 
        x2n = np.random.uniform(-1*sqrt(dispersal), sqrt(dispersal))    
        transmission_mat[meta_p][newpop] = ["{},{}".format(int(transmission_mat[meta_p][dpop][0].split(",")[0]) + x1n, int(transmission_mat[meta_p][dpop][0].split(",")[1]) + x2n),[]]
        for tpop in range(1,len(transmission_mat[meta_p].keys())):
            dist = sqrt((float(transmission_mat[meta_p]["pop_{}".format(tpop)][0].split(",")[0]) - x1n)**2 + (float(transmission_mat[meta_p]["pop_{}".format(tpop)][0].split(",")[1]) - x2n)**2)        
            if dist < dispersal:
                dist = 1
            else: 
                dist = 0
            transmission_mat[meta_p][newpop][1].append(dist)
            transmission_mat[meta_p]["pop_{}".format(tpop)][1].append(dist)
        transmission_mat[meta_p][newpop][1].append(1)   
        #if pop dies it is N in the transmission_mat, it is removed from meta_popdict
         
        return transmission_mat
    
#-----------------------------------------------------------------------------
         
class wbsims:
    def __init__(self,args):
        self.numberGens = self.numberGens
        
    def wb_cont(self,numberGens):
        '''this will call other functions to intialize, then functions from the life cycle
        num_gen: how long to run
        prevalance: .80 is 80%
        hostpop_size: number of potentially infected hosts'''
        f = open("sims.out",'w')
        f.write("time\tvillage\tnpops\tnumA1\tnumA2\tnumA2\tnumA4\tnumA5\tnumA6\tnumA7\tnumA8\tnumMF\ttrans_event\tRzero\tnum_uniqHapsAdult\tnum_uniqHapsMF\n")
        
        #initialize
        meta_popdict, hap_pop, transmission_mat, dispersal = self.worm_burden([.8], [100])  
        #set counters
        time_hours = 0 #hour counter for months
        bites_person = 0 #counter for bites
        MFpos_mosq = 0 #counter for infected mosquitoes
        month = 0 #keep track of months for maturation
        #set empty recording variables
        sum_adults = []
        #how long to run the simulation
        sim_time = numberGens * 10 * 720 #gens per year are 1.2, 12/1.2 is 10 months for 1 generation    
        
        #gillespie: rates/hour
        ##for 1 village
        events = ["bite", "infbite","L3_transmission"] #what can happen
        #rate of bite is: 20 bites per person per month * hostpopsize / 720 = bites/hour; rate_b = (20.00*hostpopsize[0])/720   
        #rate of bites that have MF: rate of bite * prev * prob pickup MF; rate_ib = rate_b * len(meta_popdict.keys()) * 0.37 
        #rate of transmission: rate_ib * number_L3 * reinjected to new host; number_L3 = 4.395(1-exp^-(.055(m))/4.395)^2; leave with bite = .414; reinjected = .32; m=(sum[donor["MF"]]/235)/50    
        #(sum[allMFstages]/235)/50 number of MF in 20ul of blood 235ml is 5% of total body blood and there are 50*20ul in 1 ML
        
        #for 2 villages there would be an extra rate category determining the rate of migration of hosts between 2 villages
        #events = ["bite1","bite2",infbite1,infbite2,L3trans1,L3trans2,Hostmig]
        while time_hours <= sim_time:
            rates = [2.7, 0.7992, 0.28] #rates of events per hour bites, Y, Z
    #        rb = 20*hostpopsize[0]/720.0
    #        rib = rb * len(meta_popdict.keys())
    #        mf = 22.5 #(sum(allMF stages)/235)/50
    #        trL3 = rib * (4.395*(1-math.exp((.055*(mf))/4.395))**2) * 0.13248
    #        rates = [rb,rib,trL3]
            random_t = np.random.uniform(0,1,len(rates)) #vector of random variables
            wait_t = -(np.log(random_t))/rates #wait time to next event
            event_t = np.argmin(wait_t) #which event has the shortest wait time
            time =+ wait_t[event_t] #update time
            if events[event_t] is "bite":
                bites_person += 1
            elif events[event_t] is "infbite":
                MFpos_mosq += 1
                bites_person += 1
            elif events[event_t] is "L3_transmission":
                meta_p = random.choice(meta_popdict.keys())
                meta_popdict, transmission_mat = self.transmission(transmission_mat,meta_p)
                bites_person += 1
                MFpos_mosq += 1
            #update times
            time_hours += time
            if time_hours >= 720*month: #each month is survival life cycle
                month += 1  #next month counter          
                meta_popdict, hap_pop, rzero_freq, hap_freq = self.maturation(meta_popdict,hap_pop,month) #update matrices    
                #sum data for month and print to out file
                for metap in meta_popdict.keys():
                    sum_adults = []                
                    sum_adults = sum_adults + [0]*9             
                    for npop in metap:
                        for i in range(0,7):
                            sum_adults[i] += len(npop["A_{}".format(i+1)])
                        for j in range(1,12):
                            sum_adults[8] += len(npop["MF_{}".format(j)])
                    f.write("%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f \t%i\n")%(month,metap,len(metap.keys()),sum_adults[0],sum_adults[1],sum_adults[2],sum_adults[3],sum_adults[4],sum_adults[5],sum_adults[6],sum_adults[7],sum_adults[8],0,np.mean(rzero_freq.values()),len(hap_freq.keys()))
            f.close()
        def wb_discrete(self,numberGens):
    
def main():
    args = get_args()
    if self.continuous: #make certain only 1 village is passed
        if len(self.hostpopsize) > 1: 
            raise ValueError("invoking continuous will only use 1 village")
        wbsims.wb_cont(self.numberGens)
    else:
        wbsims.wb_discrete(self.numberGens)
    
if __name__ == '__main__':
    main()
    
##Extra infor for worms
    #   return output_table.txt, sample.ms.out #in ggplot2 format for easy plotting, ms style with the first seq all 0's to pipe through ms2dna     
    # -t 2358.4794411347116 -eN 6.6821231193e-06 0.01572188706250533 -eN 1.42572368508e-05 0.00040268391361331134 -eN 2.2844591757e-05 0.00035830347485127096 -eN 3.2579677677e-05 0.00033732875077569337 -eN 4.3615813734e-05 0.00037284734435455506 -eN 5.6126628747000005e-05 0.00046205198060262633 -eN 7.0309283646e-05 0.0005462373570481679 -eN 8.6387439486e-05 0.0005799713845657212 -eN 0.000104614225461 0.0005540858410826578 -eN 0.000125276902926 0.0005296681469822175 -eN 0.00017525486667000002 0.0006835004731284024 -eN 0.00023948268963 0.0010190469552510712 -eN 0.00032202527898000005 0.001383597270671466 -eN 0.0004281042193500001 0.0016691809241864753 -eN 0.00056442934239 0.0020381379966736206 -eN 0.0007396248487800001 0.002753402774170087 -eN 0.0009647741508000001 0.003985880273745958 -eN 0.00125412371565 0.005853123192477116 -eN 0.00162597558966 0.00812486945657527 -eN 0.0021038555238 0.009899002633495827 -eN 0.0027180012207 0.010009088398246529 -eN 0.0035072597436 0.007449897915795466 -eN 0.004521557328 0.002036698834187558 -eN 0.0058250666766 1.0000330198024472
    '''at each year collect all the seq from the MF and store them. at end transform them back to binary '001001010'
    meta1_pop1 : 'seq_base[1]','seq_base[4]'
    meta1_pop2 : 
    
    ---if the rate here is the mean of the Poisson then this describes the expected number of events per month
    ---the waiting time to the first event is an exponential with rate and time
    ---the waiting time to the nth event is a gamma distribution
    
    -each person receives 10 bites per month from mosquitoes, proportion of mosquitoes that pick up MF when biting infected host .37; 10 bites/month per person on avg
        3.7 of those bites pick up MF; 1 in every 2.7
    -proportion of these bites that are anopheles vs culex will determine which uptake function to use
        -uptake rate in mosquito determines the number of L3 produced from number of MF in 20ul blood. Nmf/50 is number in 20ul blood
        culex: L3 = Ks1(1-exp^-(r1*m/Ks1))
            where L3 is the number of L3 produced; m is the number of MF per 20ul, r1 is the rate of development with more MF ingested
            Ks1 is the max limiting value of L3 developing.
        anopheles: L3 = Ks2(1-exp^-(r2(m-T))/Ks2)^2
            where all is the same as culex except it is squared and T is the threshold density where above this the MF get a facilitation
            effect.
        values from gambhir michael 2008: 
            Ks1: 4.406 (+- 0.362)
            r1: 0.019 (+- .058)
            Ks2: 4.395 (+- 0.332)
            r2: 0.055 (+- 0.004)
            T: 0
    proportion of L3 that leave mosquito with bite 0.414; 1 in every 2.4
    proportion of L3 that enter the host 0.32; 1 in every 3.125
    
    #unused gillespie:
        #random_t = np.random.uniform(0,1,len(rates)) #vector of random variables
        #wait_t = -(np.log(random_t))/rates #wait time to next event
        #event_t = np.argmin(wait_t) #which event has the shortest wait time
        #time =+ wait_t[event_t] #update time
        #rate of bite is: 20 bites per person per month * hostpopsize / 720 = bites/hour; rate_b = (20.00*hostpopsize[0])/720   
        #rate of bites that have MF: rate of bite * prev * prob pickup MF; rate_ib = rate_b * len(meta_popdict.keys()) * 0.37 
        #rate of transmission: rate_ib * number_L3 * reinjected to new host; number_L3 = 4.395(1-exp^-(.055(m))/4.395)^2; leave with bite = .414; reinjected = .32; m=(sum[donor["MF"]]/235)/50    
        #(sum[allMFstages]/235)/50 number of MF in 20ul of blood 235ml is 5% of total body blood and there are 50*20ul in 1 ML

    '''
    
