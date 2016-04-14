# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:56:02 2016

@author: stsmall
"""
#import all that is needed
import numpy as np
from collections import defaultdict 
from math import pi, exp, sqrt
import subprocess, re, random

#### BEGIN initialize functions
def migration_matrix(villages,distance_m):
    '''calculated a migration matrix for use with ms (hudson 2000) from euclidian
    distances. The value of 2Nm is weighted by the distance from the next population
    as an exponential random variable. The highest this can be is 2Nm/1'''    
    m = .1 #10% of migrants
    Ne = 10000
    mig = []
    for meters in distance_m:
        mig.append((m)/(np.random.exponential(meters)))
    return "0,{},{},0".format(4*Ne*mig[0],4*Ne*mig[0]) #mig_matrix is symmetrical 
            #"A,{}"
            #"{},B"
    
    
def ms_outcall(meta_popsize, theta, distance_m):
    '''external call to ms (Hudson 2000) using the migration matrix calculated
above from the euclidian distances. This function will then read in the stdout from
ms. The continious location of Seg given by ms is transformed to discrete coordinates. 
It then builds a list of lists from the discrete coordinates which is then pointer
referenced in the dictionary recording the metapopulation information
    metapop_init = [100,100] #2 metapops
    theta = [metapop1,metapop2]    
    ''' 
    if len(meta_popsize) == 1: #1 village
        mscmd = "./ms 1 {} -t {} | tee temp_file".format(meta_popsize[0],theta[0])        
    else: #set up for 2 villages
        num_subpops = len(meta_popsize) #-I num_pops
        total_inds = sum(sum(meta_popsize,[])) # how many
        sub_pop = " ".join(map(str,map(sum,meta_popsize))) #-I X i j ...
        t = 0.05 #2000gens/4*N0         
        mscmd = "./ms 1 {} -t {} -I {} {} -n 1 .5 -n 2 .5 -ma {} -ej {} 1 2| tee temp_file".format(total_inds,theta[0],num_subpops,sub_pop,migration_matrix(len(meta_popsize),distance_m),t)       
        #-ej the villages split at 2000 gens in the past which is 2000 years
        #after the split the pops are half the size of ancestral and share migrants at rate -ma
        #pop continues to go back until ancestral time    
    subprocess.Popen(mscmd, shell=True)     
    #parses ms output    
    hap_pop = []  
    with open("temp_file",'r') as ms:
        for line in ms:
            if line.startswith("positions"):
                line = ms.next()                
                for line in ms:            
                    hap = line.rstrip("\n")                    
                    hap_pop.append([m.start() for m in re.finditer("1", hap)])    
    
    return hap_pop
    
def worm_burden(prevalance,hostpop_size,distance_m):
    '''
    num_pops = [15,25] #number of infections in meta1 and meta2; prevalance * population
    mu = [5,1] #avg_burden, average number of adult female worms in infections;
    size = [1,50] #size = dispersion, size parameter for negative binomial distribution
    hap_pop is from ms_outcall    
    '''  
    # number of adult worms per infection    
    num_inf = []
    for i,j in zip(prevalance,hostpop_size):
        num_inf.append(round(i*j))
    mu = [5,1]
    sigma = [1,50]
    #worm burden    
    pop_init=[]
    for i,j,k in zip(mu,size,num_inf):
        wb_burden = np.random.negative_binomial(i,i/float(i+j),k) # number of successes, prob of success (size/size+mu),number of values to return
        pop_init.append(np.array(wb_burden).tolist())
    
    #make populations with haplotypes
    meta_popsize = []
    for inf in pop_init:
        meta_popsize.append(sum(inf))
    hap_pop = ms_outcall(meta_popsize, theta = [15.2], distance_m)  
   
   #meta_popdict[meta1][pop][age][list]
    meta_popdict = defaultdict(lambda : defaultdict(lambda : defaultdict(list)))            
    pop = 1 #initial value
    meta = 1 #initial value
    k = 0 #initial start of happop
    for metapop in pop_init:
        for wb_a1 in metapop:
            j = 0 #initializes         
            while j < wb_a1:
                meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"].append(random.uniform, k)
                j += 1 #counts the A1 in population
                k += 1 #counts the haps in hap_pop[]
            pop += 1 #advances the pop counter
        meta += 1 #advances the meta counter
        pop = 1 #resets pop counter for new meta
    return meta_popdict, hap_pop #dictionary, max(max(hap_pop))

def trans_init(prevalance,hostpop_size):
    '''initializes locations for above infections within
    prevalence = [.8,.2]
    hostpop_size = [100,100]    
    '''
    from sklearn.metrics.pairwise import euclidean_distances
    sigma = 2
    dispersal = 4*pi*sigma**2
    size = 1
    mu = 100
    metan = 1
    transmission_mat = defaultdict(lambda : defaultdict(dict))      
    for i,j in zip(prevalance, hostpop_size):    
        x1 = np.random.negative_binomial(size,size/float((size+mu)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
        x2 = np.random.negative_binomial(size,size/float((size+mu)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
        X = np.vstack((x1,x2)).T    
        dist = euclidean_distances(X, X) #this is a sparse matrix   
        for inf in range(1,len(x1)):
            transmission_mat["meta_{}".format(metan)]["pop_{}".format(inf)] = ["{},{}".format(x1[inf-1],x2[inf-1]), [1 if item < dispersal else 0 for item in dist[inf,:]]]
        metan += 1
    return transmission_mat, dispersal
    

#### END intialize functions    
    
### BEGIN life cycle functions

def maturation(meta_popdict, hap_pop, month,bp):
    '''step all individuals forward, after deaths. Kill them then move, then births'''
    #adult worms die at rate .01 per month or survive at .99
    #MF die at rate .10 per month or survive at .90
    #the prob that a juv worm survives to adult is .2, or .8177651 survival per month 

    #density dependent mortality    
#    b = 1.0/K #K is carrying capacity
#    S = 0.99
#    a = (S*b)/exp(-1) #where S is maximum survival
#    mort_A = a * sum(len(A_1 ... A_8)) * exp(-b * sum(len(A_1 ... A_8))) #Ricker fx 
    #fixed survival    
    mort_A = 0.8863848717161292 #year
    mort_J = 0.865880173963825 #month
    mort_M = 0.90 #month
    
    #density dep fecundity   
#    b = 1.0/K #K is carrying capacity
#    S = 0.99
#    a = (S*b)/exp(-1) #where S is maximum survival
#    fecund = 20 * a * sum(len(A_1 ... A_8)) * exp(-b * sum(len(A_1 ... A_8))) #Ricker fx 
    fecund = 20     

    #mutation
    bp = bp
    pmut = 7.6E-8    
    frac_gen = 0.125
    
    rzero_freq = []
    hap_freq = []
   #since this month to month 
    if month%12 is 0: #move adults 
        for mpop in meta_popdict.keys(): #villages
            for npop in mpop.keys(): #inf
                npop["A_8"] = random.sample(npop["A_7"],int(round(len(npop["A_7"])*mort_A)))
                npop["A_7"] = random.sample(npop["A_6"],int(round(len(npop["A_6"])*mort_A)))
                npop["A_6"] = random.sample(npop["A_5"],int(round(len(npop["A_5"])*mort_A)))
                npop["A_5"] = random.sample(npop["A_4"],int(round(len(npop["A_4"])*mort_A)))
                npop["A_4"] = random.sample(npop["A_3"],int(round(len(npop["A_3"])*mort_A)))
                npop["A_3"] = random.sample(npop["A_2"],int(round(len(npop["A_2"])*mort_A)))
                npop["A_2"] = random.sample(npop["A_1"],int(round(len(npop["A_1"])*mort_A)))
                # first year adults, later add from different age classes with maturity an evolvable trait            
                npop["A_1"] = random.sample(npop["J_12"],int(round(len(npop["J_12"])*mort_J)))               
                # juvenille
                npop["J_12"] = random.sample(npop["J_11"],int(round(len(npop["J_11"])*mort_J)))
                npop["J_11"] = random.sample(npop["J_10"],int(round(len(npop["J_10"])*mort_J)))
                npop["J_10"] = random.sample(npop["J_9"],int(round(len(npop["J_9"])*mort_J)))
                npop["J_9"] = random.sample(npop["J_8"],int(round(len(npop["J_8"])*mort_J)))
                npop["J_8"] = random.sample(npop["J_7"],int(round(len(npop["J_7"])*mort_J)))
                npop["J_7"] = random.sample(npop["J_6"],int(round(len(npop["J_6"])*mort_J)))
                npop["J_6"] = random.sample(npop["J_5"],int(round(len(npop["J_5"])*mort_J)))
                npop["J_5"] = random.sample(npop["J_4"],int(round(len(npop["J_4"])*mort_J)))
                npop["J_4"] = random.sample(npop["J_3"],int(round(len(npop["J_3"])*mort_J)))
                npop["J_3"] = random.sample(npop["J_2"],int(round(len(npop["J_2"])*mort_J)))
                npop["J_2"] = random.sample(npop["J_1"],int(round(len(npop["J_1"])*mort_J)))
                # microfilaria                
                npop["MF_12"] = random.sample(npop["MF_11"],int(round(len(npop["MF_11"])*mort_M)))
                npop["MF_11"] = random.sample(npop["MF_10"],int(round(len(npop["MF_10"])*mort_M)))
                npop["MF_10"] = random.sample(npop["MF_9"],int(round(len(npop["MF_9"])*mort_M)))
                npop["MF_9"] = random.sample(npop["MF_8"],int(round(len(npop["MF_8"])*mort_M)))
                npop["MF_8"] = random.sample(npop["MF_7"],int(round(len(npop["MF_7"])*mort_M)))
                npop["MF_7"] = random.sample(npop["MF_6"],int(round(len(npop["MF_6"])*mort_M)))
                npop["MF_6"] = random.sample(npop["MF_5"],int(round(len(npop["MF_5"])*mort_M)))
                npop["MF_5"] = random.sample(npop["MF_4"],int(round(len(npop["MF_4"])*mort_M)))
                npop["MF_4"] = random.sample(npop["MF_3"],int(round(len(npop["MF_3"])*mort_M)))
                npop["MF_3"] = random.sample(npop["MF_2"],int(round(len(npop["MF_2"])*mort_M)))
                npop["MF_2"] = random.sample(npop["MF_1"],int(round(len(npop["MF_1"])*mort_M)))
                #count A_1                
                rzero_freq.append([i[0]for i in npop["A_1"]])
                #reset A_1
                for subl in npop["A_1"]:
                    subl[0] = random.random()
                # MF_1
                mf = []                
                for i in range(1,8):
                    for j in npop["A_{}".format(i)]:
                        mf.append([j]*np.random.poisson(fecund)) 
                #set up mutation        
                num_muts = np.random.binomial(len(mf), bp * pmut * frac_gen)
                if num_muts != 0:
                    mf,hap_pop = mutation(mf,hap_pop)
                npop["MF1"] = mf
                hap_freq.append([i[1]for i in mf])
                
    else: #not adults
        for mpop in meta_popdict.keys():        
            for npop in mpop.keys(): #inf
                # first year adults, later add from different age classes with maturity an evolvable trait            
                npop["A_1"].append(random.sample(npop["J_12"],int(round(len(npop["J_12"])*mort_J))))
                # juvenilles                
                npop["J_12"] = random.sample(npop["J_11"],int(round(len(npop["J_11"])*mort_J)))
                npop["J_11"] = random.sample(npop["J_10"],int(round(len(npop["J_10"])*mort_J)))
                npop["J_10"] = random.sample(npop["J_9"],int(round(len(npop["J_9"])*mort_J)))
                npop["J_9"] = random.sample(npop["J_8"],int(round(len(npop["J_8"])*mort_J)))
                npop["J_8"] = random.sample(npop["J_7"],int(round(len(npop["J_7"])*mort_J)))
                npop["J_7"] = random.sample(npop["J_6"],int(round(len(npop["J_6"])*mort_J)))
                npop["J_6"] = random.sample(npop["J_5"],int(round(len(npop["J_5"])*mort_J)))
                npop["J_5"] = random.sample(npop["J_4"],int(round(len(npop["J_4"])*mort_J)))
                npop["J_4"] = random.sample(npop["J_3"],int(round(len(npop["J_3"])*mort_J)))
                npop["J_3"] = random.sample(npop["J_2"],int(round(len(npop["J_2"])*mort_J)))
                npop["J_2"] = random.sample(npop["J_1"],int(round(len(npop["J_1"])*mort_J)))
                # microfilaria                
                npop["MF_12"] = random.sample(npop["MF_11"],int(round(len(npop["MF_11"])*mort_M)))
                npop["MF_11"] = random.sample(npop["MF_10"],int(round(len(npop["MF_10"])*mort_M)))
                npop["MF_10"] = random.sample(npop["MF_9"],int(round(len(npop["MF_9"])*mort_M)))
                npop["MF_9"] = random.sample(npop["MF_8"],int(round(len(npop["MF_8"])*mort_M)))
                npop["MF_8"] = random.sample(npop["MF_7"],int(round(len(npop["MF_7"])*mort_M)))
                npop["MF_7"] = random.sample(npop["MF_6"],int(round(len(npop["MF_6"])*mort_M)))
                npop["MF_6"] = random.sample(npop["MF_5"],int(round(len(npop["MF_5"])*mort_M)))
                npop["MF_5"] = random.sample(npop["MF_4"],int(round(len(npop["MF_4"])*mort_M)))
                npop["MF_4"] = random.sample(npop["MF_3"],int(round(len(npop["MF_3"])*mort_M)))
                npop["MF_3"] = random.sample(npop["MF_2"],int(round(len(npop["MF_2"])*mort_M)))
                npop["MF_2"] = random.sample(npop["MF_1"],int(round(len(npop["MF_1"])*mort_M)))
                #count A_1                
                rzero_freq.append([i[0]for i in npop["A_1"]])
                #reset A_1
                for subl in npop["A_1"]:
                    subl[0] = random.random()
                # MF_1
                mf = []                
                for i in range(1,8):
                    for j in npop["A_{}".format(i)]:
                        mf.append([j]*np.random.poisson(fecund)) 
                #set up mutation        
                num_muts = np.random.binomial(len(mf), bp * pmut * frac_gen)
                if num_muts != 0:
                    mf,hap_pop = mutation(mf,hap_pop)
                npop["MF1"] = mf
                hap_freq.append([i[1]for i in mf])
                
    return meta_popdict, hap_pop, {x:rzero_freq.count(x) for x in rzero_freq},{y:hap_freq.count(y) for y in hap_freq}
    
def mutation(mf,hap_pop):
   '''this is run every time the prob of a mutation is true, updates seq_base'''   
   #calculate the number of mutations expected
   mut_mf = [random.randrange(len(mf)) for i in range(num_muts)] #choose random index in mf list
       #add new sequence to hap_pop
   for m in mut_mf:                 
       new_hap = copy.copy(hap_pop[mf[m][1]])
       new_hap.append(max([max(a) for a in zip(hap_pop)]))
       hap_pop.append(new_hap)
       mf[m][1] = len(hap_pop)-1
   return mf, hap_pop

def transmission(transmission_mat,metap,dispersal):    
    '''this has a continious prob using a gillepsie algorithm to get the waiting time between transmission events'''    
    #this executes if tranmission is true     
    #pick a random donating pop:
    dpop = random.choice(meta_popdict[meta_p].keys())    
    r = random.uniform(0,1)    
    #new infection    
    if r > 1 - (1/transmission_mat[metap][dpop][1]):
        npop = "pop_{}".format(len(meta_popdict[metap].keys()))
        meta_popdict[metap][npop]["J_1"].append(random.choice(meta_popdict[metap][dpop]["MF_{}".format(random.uniform(1,12))]))
        transmission_mat = new_infection(metap,dpop,npop,dispersal) #find new transmission position
    #reinfection
    else:     
        rpop = random.choice([i for i, x in enumerate(transmission_mat[metap][dpop][1]) if x == 1]) #choose a random value that is 1 as the receiving pop        
        meta_popdict[metap]["pop_{}".format(rpop)]["J_1"].append(random.choice(meta_popdict[metap][dpop]["MF_{}".format(random.uniform(1,12))]))
    return meta_popdict, transmission_mat

def new_infection(metap,dpop,npop,dispersal):
    '''this is run every time a new individual is infected to rebuild the transmission prob matrix, updates'''        
    x1n = np.random.uniform(-1*sqrt(dispersal), sqrt(dispersal)) 
    x2n = np.random.uniform(-1*sqrt(dispersal), sqrt(dispersal))    
    transmission_mat[metap][npop] = ["{},{}".format(int(transmission_mat[metap][dpop][0].split(",")[0]) + x1n, int(transmission_mat[metap][dpop][0].split(",")[1]) + x2n),[]]
    for tpop in range(1,len(transmission_mat[metap].keys())):
        dist = sqrt((int(transmission_mat[metap]["pop_{}".format(tpop)][0].split(",")[0]) - x1n)**2 + (int(transmission_mat["meta_1"]["pop_" + tpop][0].split(",")[1]) - x2n)**2)        
        if dist > dispersal:
            dist = 1
        else: 
            dist = 0
        transmission_mat[metap][npop][1].append(dist)
        transmission_mat[metap]["pop_{}".format(tpop)][1].append(dist)
    #if pop dies it is N in the transmission_mat, it is removed from meta_popdict
     
    return transmission_mat

##### END life cycle functions    
   
##### BEGIN sims   
    
def wb_sims(num_gens, prevalance, hostpop_size, distance_m,bp):
    '''this will call other functions to intialize, then functions from the life cycle
    num_gen: how long to run
    prevalance: .80 is 80%
    hostpop_size: number of potentially infected hosts'''
    f = open("sims.out",'w')
    f.write("time\tvillage\tnpops\tnumA1\tnumA2\tnumA2\tnumA4\tnumA5\tnumA6\tnumA7\tnumA8\tnumMF\ttrans_event\tRzero\tnum_uniqHapsAdult\tnum_uniqHapsMF\n")
    #initialize
    meta_popdict, hap_pop = worm_burden(prev,hostpop_size,distance_m) #[.8,.1] [100,1000] [2000]   
    transmission_mat, dispersal = trans_init(prev, hostpop_size)
    gens = 0
    gen = 0
    month = 1
    sum_adults = []    
    #iterate
    #gillespie
    #rates/hour
    rates = [2.5,2.6,.001]    
    wait_t = np.random.uniform(0,1,len(rates))
    event = np.argmin(rates*wait_t)
    
    #bites/person/hour 2000/month 720 hours/month = 1.25 bites per hour
    #bites per hour
    #bites on infet = prev 
    #prob of MF = number of MF in inf
    #prob MF-L3 = number of MF in inf * maturity

    while gens <= num_gens*720:    
        wait_time = -math.log(random.uniform(0,1)/1)       
        if random.uniform(0,1) <= bite_trans:
            metap = random.choice(["meta_1","meta_2"])
            meta_popdict, transmission_mat = transmission(transmission_mat,metap,dispersal)
        else:
            bites_person += 1
        gens += wait_time
        gen += wait_time        
        if gen >= 720*month: #hours per month
            month += 1            
            meta_popdict, hap_pop, rzero_freq, hap_freq = maturation(meta_popdict,hap_pop,month,bp)    
            for metap in meta_popdict.keys():
                sum_adults = []                
                sum_adults = sum_adults + [0]*9             
                for npop in metap:
                    for i in range(0,7):
                        sum_adults[i] += len(npop["A_{}".format(i+1)])
                    for j in range(1,12):
                        sum_adults[8] += len(npop["MF_{}".format(j)])
            f.write("%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f \t%i\n")%(month,metap,len(metap.keys()),sum_adults[0],sum_adults[1],sum_adults[2],sum_adults[3],sum_adults[4],sum_adults[5],sum_adults[6],sum_adults[7],sum_adults[8],0,np.mean(rzero.values()),len(hap_freq.keys()))
    f.close()

#   return output_table.txt, sample.ms.out #in ggplot2 format for easy plotting, ms style with the first seq all 0's to pipe through ms2dna     
# -t 2358.4794411347116 -eN 6.6821231193e-06 0.01572188706250533 -eN 1.42572368508e-05 0.00040268391361331134 -eN 2.2844591757e-05 0.00035830347485127096 -eN 3.2579677677e-05 0.00033732875077569337 -eN 4.3615813734e-05 0.00037284734435455506 -eN 5.6126628747000005e-05 0.00046205198060262633 -eN 7.0309283646e-05 0.0005462373570481679 -eN 8.6387439486e-05 0.0005799713845657212 -eN 0.000104614225461 0.0005540858410826578 -eN 0.000125276902926 0.0005296681469822175 -eN 0.00017525486667000002 0.0006835004731284024 -eN 0.00023948268963 0.0010190469552510712 -eN 0.00032202527898000005 0.001383597270671466 -eN 0.0004281042193500001 0.0016691809241864753 -eN 0.00056442934239 0.0020381379966736206 -eN 0.0007396248487800001 0.002753402774170087 -eN 0.0009647741508000001 0.003985880273745958 -eN 0.00125412371565 0.005853123192477116 -eN 0.00162597558966 0.00812486945657527 -eN 0.0021038555238 0.009899002633495827 -eN 0.0027180012207 0.010009088398246529 -eN 0.0035072597436 0.007449897915795466 -eN 0.004521557328 0.002036698834187558 -eN 0.0058250666766 1.0000330198024472
'''at each year collect all the seq from the MF and store them. at end transform them back to binary '001001010'
meta1_pop1 : 'seq_base[1]','seq_base[4]'
meta1_pop2 : 

---if the rate here is the mean of the Poisson then this describes the expected number of events per month
---the waiting time to the first event is an exponential with rate and time
---the waiting time to the nth event is a gamma distribution

each population receives 10 bites per month from mosquitoes
    proportion of these are anopheles vs culex will determine which uptake fx to use
proportion of mosquitoes that pick up MF when biting infected host .37; 10 bites/month per person on avg
    3.7 of those bites pick up MF; 1 in every 2.7
uptake rate in mosquito determines the number of L3 produced from number of MF in 20ul blood. Nmf/50 is number in 20ul blood
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

'''