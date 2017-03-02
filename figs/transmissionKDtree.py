#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import ipdb
import math
import random
import pickle

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .agehost import agehost_fx

deathdict = pickle.load(open('../figs/data/acttable.p', "rb"))

def vectorbite_fx(vill,
                  month,
                  village,
                  densitydep_uptake,
                  avgMF):

    '''Counts the number of successful infectious mosquito bites

    Parameters
    --------
    bitesPperson : int
        rate of bites per person per unit time, here hours
    hours2bite : int
        number of hours mosq bite per night/day
    hostpopsize : int
        total population of village
    prev_t : float
        percent of people infected
    densitydep_uptake : Boolean
        use the density dependece function for developing L3
    avgMF : float
         average number of MF per host per village
    bednets : Boolean, list
        use values from bednets
    bnstart : int, list
        month to start bednets
    bnstop : int, list
        mont to stop bednets
    bncoverage : int, list
        percent of population using bednets
    month : int
        current time, month
    Returns
    ------
    L3trans : int
         all new juv are age class 0
    '''

    bitesPperson = village[vill].bpp
    hours2bite = village[vill].h2b
    hostpopsize = village[vill].hostpopsize
    prev_t = village[vill].prev
    bednets = village[vill].bn
    bnstart = village[vill].bnstr
    bnstop = village[vill].bnstp
    bncoverage =  village[vill].bncov
    #print("vectorbite")
    if bednets:
        if month > bnstart and month < bnstop:
             totalbites = ((1 - bncoverage) * bitesPperson * hours2bite * 30
                                   * hostpopsize)
        else:
             totalbites = (bitesPperson * hours2bite * 30
                                   * hostpopsize)
    else:
        totalbites = (bitesPperson * hours2bite * 30
                              * hostpopsize)
    # 0.37 is prob of bite on infected host picking up MF
    infbites = np.random.binomial(totalbites, (prev_t * 0.37))
    if densitydep_uptake: #values for anopheles from CITE
       #number of MF in 20ul of blood
       #235ml is 5% of total host blood
       #there are 50 units of 20ul in 1ml
       mfBlood = avgMF / 50.0
       # 0.414 is proportion of  L3 that leave mosquito per bite
       # 0.32 proportion of L3 that make it into the host
       L3trans = round(infbites * (4.395 * (1 - math.exp( -(0.055 * (mfBlood))
           / 4.395)) ** 2) * (0.414 * 0.32))
    else:
       # 0.414 is proportion of  L3 that leave mosquito per bite
       # 0.32 proportion of L3 that make it into the host
       L3trans = np.random.binomial(infbites, (0.414 * 0.32))
    return(int(L3trans))


def new_infection_fx(dispersal,
                     mfhostidx,
                     dfHost):
    '''A transmission event infecting a naive host

    Parameters
    ----------
    dispersal: float
         dispersal distance as 2*sigma
    transMF
        row of transmitted MF
    dfHost: df
        host dataframe

    Returns
    ------
    dfHost: df
         updated with new host information
    newhostidx: str
          new host index
    '''
    print("newfx")
    #how close
    mindist = 5
    #how far
    maxdist = dispersal
    #assign new position
    maxc = int( math.sqrt(( maxdist ** 2 / 2 )))
    #makes the below statement true
    newpts = dfHost[dfHost.hostidx == mfhostidx].coordinates.values[0]

    while next((True for elem in dfHost['coordinates'].values if np.array_equal(elem, newpts)),False) is True:
         #shift of new points
         x_new = random.randint( -maxc, maxc )
         y_new = random.randint( -maxc, maxc )
         #test that new shifts are min distance
         if math.sqrt(x_new ** 2 + y_new ** 2) > mindist:
              newhostptX = dfHost.coordinates[dfHost.hostidx == mfhostidx].values[0][0] + x_new
              newhostptY = dfHost.coordinates[dfHost.hostidx == mfhostidx].values[0][1] + y_new
              newpts = np.array([newhostptX, newhostptY])
    #copy village
    vill = mfhostidx[:mfhostidx.rfind('h')][-1]
    #new host index
    old_hostidx = dfHost[dfHost.village == vill].hostidx.iloc[-1]
    new_hostidx = old_hostidx[:old_hostidx.rfind('h')] + 'h' + str(int(old_hostidx.split('h')[1]) + 1)
    #radnom sex
    sex = random.choice("01")
    #age and agedeath function from wbsims_initialize
    age, agedeath = agehost_fx(sex, deathdict)
    #add to dfHost at bottom
    newhostlist = [[vill, new_hostidx, sex, age, agedeath, newpts, 0, 0]]
    dfHost = pd.concat([dfHost, pd.DataFrame(newhostlist,columns=dfHost.columns)],ignore_index=True)
    return(dfHost, new_hostidx)

def transmission_fx(month,
                    village,
                    sigma,
                    densitydep_uptake,
                    dfHost,
                    dfworm):
    '''Transmission events resolved as either reinfection or new infection

    Parameters
    --------
    villages : int
        number of villages
    hostpopsize : int, list
        total population size of village
    sigma : float
        dispersal mean
    bitesPperson : int, list
        rate of bites per person per unit time, here hours
    hours2bite : int, list
        number of hours mosq bite per night/day
    densitydep_uptake : Boolean
        use the density dependece function for developing L3
    bednets : Boolean, list
        use values from bednets
    bnstart : int, list
        month to start bednets
    bnstop : int, list
        mont to stop bednets
    bncoverage : int, list
        percent of population using bednets
    month : int
        current time, month
    dfHost: df
         chooses the donating and accepting hosts for transmission
    dfMF : df
         donating MF genotype to transmit
    dfJuv : df
         donated MF move to Juvenille age class

    Returns
    ------
    dfJuv : df
        all new juv are age class 0
    dfMF : df
        removes MF that were transmitted
    dfHost : df
        in the case of newinfection, new host
    L3trans : int
        number of transmitted MF
    '''
    print("transmission_fx")

    mfiix = dfworm.meta[dfworm.meta.stage == "M"].index.values
    if dfworm.meta.ix[mfiix].shape[0] > 0:
        assert dfworm.pos['0'].shape[0] == dfworm.h1['0'].shape[1]
    else: pass
    dispersal = 2 * sigma
    new_hostidx = []
    new_juv = []
    tree = cKDTree(np.vstack(dfHost.coordinates))
    mfiix_vill = np.array([dfworm.meta.ix[mfiix][dfworm.meta.village == vill].index.values for vill in range(len(village))])
    for vill in range(len(village)):
        infhost = (dfHost.village == vill).sum()
        prev_t = infhost / float(village[vill].hostpopsize)
        village[vill].prev = prev_t
        avgMF = mfiix_vill[vill].shape[0]/float(infhost)
        L3trans = vectorbite_fx(vill, month, village, densitydep_uptake, avgMF)
        print("village is %i transmitted is %i" %(vill, L3trans))
        if L3trans != 0:
            if L3trans > mfiix_vill[vill].shape[0]:  #more transmision events than possible MF
                transMF = mfiix_vill[vill]
            else:
                transMF = np.random.choice(mfiix_vill[vill], L3trans, replace=False)
            tcount = ''
            new_juv.extend(transMF)
            for idx in transMF:
                if dfworm.meta.ix[idx].hostidx != tcount:
                    mfhostidx = dfworm.meta.ix[idx].hostidx
                    transhostidx = dfHost[dfHost.hostidx == mfhostidx].index[0] #index of donating host
                    transhost = tree.query_ball_point(dfHost.ix[transhostidx].coordinates, dispersal)
                    tcount = mfhostidx
                else:
                    pass
                if len(dfHost[dfHost.village == vill]) < village[vill].hostpopsize:
                     prob_newinfection = 1.0 / (len(transhost) + 1)
                else: #everyone is already infected
                     prob_newinfection = 0
                if np.random.random() < prob_newinfection:
                    dfHost, rehostidx = new_infection_fx(dispersal, mfhostidx, dfHost)
                    new_hostidx.append(rehostidx)
                    #new host so have to resort and rebuild KDTree
                    tree = cKDTree(np.vstack(dfHost.coordinates))
                    tcount = ''
                else:
                    try: #allow self infection
                        rehostidx = transhost[np.random.randint(len(transhost) + 1)]
                    except IndexError:
                        rehostidx = transhostidx
                    new_hostidx.append(dfHost.ix[rehostidx,'hostidx'])

        else:
            print("dfMF is empty")
    dfworm.meta.ix[new_juv].stage = "J"
    dfworm.meta.ix[new_juv].hostidx = new_hostidx
#    ipdb.set_trace()
    return(village, dfHost, dfworm, L3trans)
