#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""

import numpy as np
from scipy.stats import weibull_min
from figs.hostmda import hostmda_fx
from figs.fecundity_mda import fecunditymda_fx
#from figs.fecundity_mda import fecunditymda_sel1_fx
#from figs.fecundity_mda import fecunditymda_sel2_fx
from figs.host_migration import hostmigration_fx

def kill_adults(dfworm, dfHost, month, shapeAdult, scaleAdult,
        village):
    """
    """
    adiix = dfworm.meta.index[dfworm.meta.stage == "A"].values
    #Adult survival is based on weibull cdf
    kill_adult_rand = np.random.random(adiix.shape[0])
    try:
        kill_adult_fx_age = weibull_min.cdf(dfworm.meta.age[adiix], shapeAdult,
                loc=0, scale=scaleAdult)
    except TypeError:
        kill_adult_fx_age = weibull_min.cdf(0, shapeAdult, loc=0, scale=scaleAdult)
    dieAdult = adiix[np.where(kill_adult_rand < kill_adult_fx_age)]
    dfworm.meta.ix[adiix, "age"] += 1 #2 - 21
    ##host survival is from act table
    dfHost = dfHost.query("age < agedeath")
    diehost = dfHost.hostidx.values

    dead_worms = np.append(dieAdult,
            dfworm.meta[~dfworm.meta.hostidx.isin(diehost)].index.values)
    dfworm.drop_worms(dead_worms)
    return(dfHost, dfworm)

def survivalmda_fx(month,
                   village,
                   surv_Juv,
                   shapeMF,
                   scaleMF,
                   shapeAdult,
                   scaleAdult,
                   fecund,
                   locus,
                   mutation_rate,
                   recombination_rate,
                   basepairs,
                   selection,
                   hostmigrate,
                   mdalist,
                   densitydep_surv,
                   densitydep_fec,
                   dfHost,
                   dfworm):

    '''base survival function
    Parameters
    ---------
    month: int
         time in months of the simulation
    surv_Juv: float
         survival prob of juvenille life stage
    shapeMF: float
         shape parameter for weibull distribution of MF
    scaleMF: int
         scale parameter for weibull distribution of MF
    shapeAdult: float
         shape parameter for weibull distribution of Adults
    scaleAdult: int
         scale parameter for weibull distribution of MF
    dfMF : df
        dataframe of MF
    dfAdult : df
        dataframe of adult worms
    dfJuv : df
        dataframe of juvenille worms
    dfHost : df
        dataframe of hosts
    basepairs : int, list
        length of loci
    hostmigrate : float
        rate of migration per year between villages
    selection : boolean
        T/F for selection
    dfSel : df
        dataframe of cds positions and fitness
    cds_coordinates : list
        list of coding seq coordinates
    mda : boolean
         mda in village or not
    mda_start : list, int
         time in months that mda began in each village
    mda_num : list, int
         how many mdas to simulate
    mda_freq : int
         how often
    mda_coverage : float
         what percentage of the population get mda
    mda_macro: float
         percent of adults killed by drug
    mda_micro: float
         percent of MF killed by drug
    mda_juvcide: float
         percent of juveniile killed by drug; could be 0 or avg micro/macro
    mda_sterile : float
         percent of adult worms perm sterilized
    mda_clear : int
         clearance time

    Returns
    -------
    dfMF
    dfAdult
    dfJuv
    dfHost
    dfSel

   '''
    mda_start = mdalist[0]
    mda_num = mdalist[1]
    mda_freq = mdalist[2]
    mda_coverage = mdalist[3]
    mda_macro = mdalist[4]
    mda_juvicide = mdalist[5]
    mda_micro = mdalist[6]
    mda_sterile = mdalist[7]
    mda_clear = mdalist[8]

    ##calculates time since mda as clear_count
    if month < mda_start:
            clear_count = 0
    elif month >= mda_start and month < (mda_start + mda_freq * mda_num + mda_freq):
        if month%mda_freq == 0:
             clear_count = 1
        elif (month - mda_clear)%mda_freq == 0:
             clear_count = mda_clear + 1
        else:
             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
                  clear_count = ((month-mda_freq)%mda_freq + 1)
             else:
                  clear_count = 0
    elif month > (mda_start + mda_freq * mda_num):
        clear_count = 0
    else:pass
    hostmda = ''
#######################
    ##assign MDA to hosts
    if clear_count == 1:
        ##add mda to dfHost per coverage
        dfHost = hostmda_fx(village, dfHost, mda_coverage)
        ##apply MDA effect to host populations of worms
        hostmda = dfHost[dfHost.MDAstate == 1].hostidx.values
        wormmda = dfworm.meta[dfworm.meta.hostidx.isin(hostmda)].index.values
        MFmda = dfworm.meta.ix[wormmda].query("stage == 'M'").index.values
        Juvmda = dfworm.meta.ix[wormmda].query("stage == 'J'").index.values
        Adultmda = dfworm.meta.ix[wormmda].query("stage == 'A'").index.values

        mfdiemda = dfworm.meta.ix[MFmda].groupby("hostidx").apply(lambda y: y.sample(frac = mda_micro).index.values)
        juvdiemda = dfworm.meta.ix[Juvmda].groupby("hostidx").apply(lambda y: y.sample(frac = mda_juvicide).index.values)
        adultdiemda = dfworm.meta.ix[Adultmda].groupby("hostidx").apply(lambda y: y.sample(frac = mda_macro).index.values)
        diemda = [i for i in mfdiemda.values] + [i for i in juvdiemda.values] + [i for i in adultdiemda.values]
        dfworm.drop_worms(np.concatenate(diemda).ravel())
    else: pass
#######################
    if month%12 == 0:
        dfHost, dfworm = kill_adults(dfworm, dfHost, month, shapeAdult, scaleAdult,
                village)

        dfHost.age = dfHost.age + 1
        if hostmigrate != 0:
            dfHost = hostmigration_fx(village, dfHost, hostmigrate)
    else: pass

    ##Juv is exponential 0.866; surv_Juv
    juviix = dfworm.meta[dfworm.meta.stage == "J"].index.values
    kill_juvrand = np.random.random(juviix.shape[0])
    dieJuv = juviix[np.where(kill_juvrand > surv_Juv)]
    dfworm.meta.ix[juviix,'age'] += 1

    ##MF is weibull cdf
    mfiix = dfworm.meta[dfworm.meta.stage == "M"].index.values
    kill_mfrand = np.random.random(mfiix.shape[0])
    try:
        kill_mffxage = weibull_min.cdf(dfworm.meta.ix[mfiix].age,shapeMF,loc=0,scale=scaleMF)
    except TypeError:
        kill_mffxage = weibull_min.cdf(0,shapeMF,loc=0,scale=scaleMF)
    dieMF = mfiix[np.where(kill_mfrand < kill_mffxage)]
    dfworm.meta.ix[mfiix, 'age'] += 1

    ##move Juv age 13 to adult age 1
#    ipdb.set_trace()
    juviix12 = dfworm.meta.ix[juviix].query('age > 12').index.values
    if any(juviix12):
        #reset age to adult
        dfworm.meta.ix[juviix12,'age'] = 1
        #increase R0net for next gen
        dfworm.meta.ix[juviix12,'R0net'] += 1
        dfworm.meta.ix[juviix12,'stage'] = "A"
    else:pass
    dfworm.drop_worms(np.append(dieJuv, dieMF))
    #fecundity calls mutation/recombination
    dfAdult_mf, dfworm = fecunditymda_fx(fecund, dfworm, locus, mutation_rate,
                                         recombination_rate, basepairs, selection,
                                         densitydep_fec, mda_sterile, clear_count, mda_clear,
                                         hostmda)
    dfAdult_mf.meta.sex = np.random.choice(['M', 'F'] , size=len(dfAdult_mf.meta))
    dfAdult_mf.meta.age = 1
    dfworm.add_worms(dfAdult_mf, dfAdult_mf.meta.index.values)

    return(dfHost, dfworm)
#
#def survivalmda_sel1_fx(month,
#                        village,
#                        surv_Juv,
#                        shapeMF,
#                        scaleMF,
#                        shapeAdult,
#                        scaleAdult,
#                        fecund,
#                        locus,
#                        mutation_rate,
#                        recombination_rate,
#                        basepairs,
#                        selection,
#                        hostmigrate,
#                        mdalist,
#                        densitydep_surv,
#                        densitydep_fec,
#                        dfHost,
#                        dfAdult,
#                        dfJuv,
#                        dfMF):
#
#
#    '''base survival function
#    Parameters
#    ---------
#    month: int
#         time in months of the simulation
#    surv_Juv: float
#         survival prob of juvenille life stage
#    shapeMF: float
#         shape parameter for weibull distribution of MF
#    scaleMF: int
#         scale parameter for weibull distribution of MF
#    shapeAdult: float
#         shape parameter for weibull distribution of Adults
#    scaleAdult: int
#         scale parameter for weibull distribution of MF
#    dfMF : df
#        dataframe of MF
#    dfAdult : df
#        dataframe of adult worms
#    dfJuv : df
#        dataframe of juvenille worms
#    dfHost : df
#        dataframe of hosts
#    basepairs : int, list
#        length of loci
#    hostmigrate : float
#        rate of migration per year between villages
#    selection : boolean
#        T/F for selection
#    dfSel : df
#        dataframe of cds positions and fitness
#    cds_coordinates : list
#        list of coding seq coordinates
#    mda : boolean
#         mda in village or not
#    mda_start : list, int
#         time in months that mda began in each village
#    mda_num : list, int
#         how many mdas to simulate
#    mda_freq : int
#         how often
#    mda_coverage : float
#         what percentage of the population get mda
#    mda_macro: float
#         percent of adults killed by drug
#    mda_micro: float
#         percent of MF killed by drug
#    mda_juvcide: float
#         percent of juveniile killed by drug; could be 0 or avg micro/macro
#    mda_sterile : float
#         percent of adult worms perm sterilized
#    mda_clear : int
#         clearance time
#
#    Returns
#    -------
#    dfMF
#    dfAdult
#    dfJuv
#    dfHost
#    dfSel
#
#   '''
#    mda_start = mdalist[0]
#    mda_num = mdalist[1]
#    mda_freq = mdalist[2]
#    mda_coverage = mdalist[3]
#    mda_macro = mdalist[4]
#    mda_juvicide = mdalist[5]
#    mda_micro = mdalist[6]
#    mda_sterile = mdalist[7]
#    mda_clear = mdalist[8]
#
#    ##calculates time since mda as clear_count
#    if month < mda_start:
#            clear_count = 0
#    elif month >= mda_start and month < (mda_start + mda_freq * mda_num + mda_freq):
#        if month%mda_freq == 0:
#             clear_count = 1
#        elif (month - mda_clear)%mda_freq == 0:
#             clear_count = mda_clear + 1
#        else:
#             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
#                  clear_count = ((month-mda_freq)%mda_freq + 1)
#             else:
#                  clear_count = 0
#    elif month > (mda_start + mda_freq * mda_num):
#        clear_count = 0
##################################
#    ##assign MDA to hosts
#    if clear_count == 1:
#        ##add mda to dfHost per coverage
#        dfHost = hostmda_fx(village, dfHost, mda_coverage)
#         ##apply MDA effect to host populations of worms
#        for index, row in dfHost[dfHost.MDA == 1].iterrows():
#            #kill MF
#            try:
#                mdarandmf = np.random.random(len(dfMF.meta.hostidx == row.hostidx))
#                mdakillmf = np.where(mdarandmf < mda_micro ** dfMF.meta[dfMF.meta.hostidx == row.hostidx]["selS"])
#                dfMF.drop_worms(mdakillmf[0])
#            except ValueError:
#                pass
#             #kill Juv
#            try:
#                mdarandjuv = np.random.random(len(dfJuv.meta.hostidx == row.hostidx))
#                mdakilljuv = np.where(mdarandjuv < mda_juvicide ** dfJuv.meta[dfJuv.meta.hostidx == row.hostidx]["selS"])
#                dfJuv.drop_worms(mdakilljuv[0])
#            except ValueError:
#                pass
#            #kill Adult
#            try:
#                mdarandad = np.random.random(len(dfAdult.meta.hostidx == row.hostidx))
#                mdakillad = np.where(mdarandad < mda_macro ** dfAdult.meta[dfAdult.meta.hostidx == row.hostidx]["selS"])
#                dfAdult.drop_worms(mdakillad[0])
#            except ValueError:
#                pass
#########################################
#    print('survival pos shape')
#    print(dfAdult.pos['1'].shape[0])
#    print(dfAdult.h1['1'].shape[1])
#    #adult worms and hosts are only evaluated per year
#    if month%12 == 0:
#        #Adult survival is based on weibull cdf
#        kill_adult_rand = np.random.random(dfAdult.meta.shape[0])
#        try:
#            kill_adultfxage = weibull_min.cdf(dfAdult.meta.age, shapeAdult,
#                    loc=0, scale=scaleAdult)
#        except TypeError:
#            kill_adultfxage = weibull_min.cdf(0, shapeAdult, loc=0, scale=scaleAdult)
#        dieAdult = np.where(kill_adult_rand < kill_adultfxage)
#        dfAdult.drop_worms(dieAdult[0])
#
#        dfAdult.meta.age = dfAdult.meta.age + 1 #2 - 21
#        ##host survival is from act table
#        dfHost = dfHost[dfHost.age < dfHost.agedeath]
#        #remove all worms with dead host.hostidx from all dataframes
#        dfAdult.drop_worms(dfAdult.meta[~dfAdult.meta.hostidx.isin(dfHost.hostidx)].index.values)
#        dfJuv.drop_worms(dfJuv.meta[~dfJuv.meta.hostidx.isin(dfHost.hostidx)].index.values)
#        dfMF.drop_worms(dfMF.meta[~dfMF.meta.hostidx.isin(dfHost.hostidx)].index.values)
#        #add 1 year to all ages of hosts
#        dfHost.age = dfHost.age + 1
#        if hostmigrate != 0:
#            dfHost = hostmigration_fx(village, dfHost, hostmigrate)
#
#    print('Post killing')
#    print(dfAdult.pos['1'].shape[0])
#    print(dfAdult.h1['1'].shape[1])
#    print('Nworms :{0!s}'.format(dfAdult.h1['1'].shape[0]))
#
#    ##Juv is exponential 0.866; surv_Juv
#    #dont include age 0 which just moved from transmission fx
#    dfJuv.meta.age += 1
#    kill_juvrand = np.random.random(dfJuv.meta.shape[0])
#    dieJuv = np.where(kill_juvrand > surv_Juv)
#    dfJuv.drop_worms(dieJuv[0])
#
#    ##MF is weibull cdf
#    kill_mfrand = np.random.random(dfMF.meta.shape[0])
#    try:
#        kill_mffxage = weibull_min.cdf(dfMF.meta.age,shapeMF,loc=0,scale=scaleMF)
#    except TypeError:
#        kill_mffxage = weibull_min.cdf(0,shapeMF,loc=0,scale=scaleMF)
#    dieMF = np.where(kill_mfrand < kill_mffxage)
#    dfMF.drop_worms(dieMF[0])
#    dfMF.meta.age = dfMF.meta.age + 1 #2 - 12
#    dfMF.drop_worms(dfMF.meta.ix[dfMF.meta.age > 12].index.values) #hard cutoff at 12 months
#
#    ##move Juv age 13 to adult age 1
#    juv_rows = dfJuv.meta[dfJuv.meta.age > 12].index.values
#    if any(juv_rows):
#        #reset age to adult
#        dfJuv.meta.ix[juv_rows,"age"] = 1
#        #increase R0net for next gen
#        dfJuv.meta.ix[juv_rows, "R0net"] += 1
#        dfAdult.add_worms(dfJuv, juv_rows)
#        dfJuv.drop_worms(juv_rows)
#
#    #fecundity calls mutation/recombination
#    dfAdult_mf, dfAdult = fecunditymda_sel1_fx(fecund, dfAdult, locus, mutation_rate,
#                                         recombination_rate, basepairs, selection,
#                                         densitydep_fec, mda_sterile, clear_count, mda_clear, dfHost)
#    dfAdult_mf.meta.sex = [random.choice("MF") for i in range(len(dfAdult_mf.meta))]
#    dfAdult_mf.meta.age = 1
#    dfMF.add_worms(dfAdult_mf, dfAdult_mf.meta.index.values)
#    return(dfHost, dfAdult, dfJuv, dfMF)
#
#def survivalmda_sel2_fx(month,
#                        village,
#                        surv_Juv,
#                        shapeMF,
#                        scaleMF,
#                        shapeAdult,
#                        scaleAdult,
#                        fecund,
#                        locus,
#                        mutation_rate,
#                        recombination_rate,
#                        basepairs,
#                        selection,
#                        hostmigrate,
#                        mdalist,
#                        densitydep_surv,
#                        densitydep_fec,
#                        dfHost,
#                        dfAdult,
#                        dfJuv,
#                        dfMF):
#
#
#    '''base survival function
#    Parameters
#    ---------
#    month: int
#         time in months of the simulation
#    surv_Juv: float
#         survival prob of juvenille life stage
#    shapeMF: float
#         shape parameter for weibull distribution of MF
#    scaleMF: int
#         scale parameter for weibull distribution of MF
#    shapeAdult: float
#         shape parameter for weibull distribution of Adults
#    scaleAdult: int
#         scale parameter for weibull distribution of MF
#    dfMF : df
#        dataframe of MF
#    dfAdult : df
#        dataframe of adult worms
#    dfJuv : df
#        dataframe of juvenille worms
#    dfHost : df
#        dataframe of hosts
#    basepairs : int, list
#        length of loci
#    hostmigrate : float
#        rate of migration per year between villages
#    selection : boolean
#        T/F for selection
#    dfSel : df
#        dataframe of cds positions and fitness
#    cds_coordinates : list
#        list of coding seq coordinates
#    mda : boolean
#         mda in village or not
#    mda_start : list, int
#         time in months that mda began in each village
#    mda_num : list, int
#         how many mdas to simulate
#    mda_freq : int
#         how often
#    mda_coverage : float
#         what percentage of the population get mda
#    mda_macro: float
#         percent of adults killed by drug
#    mda_micro: float
#         percent of MF killed by drug
#    mda_juvcide: float
#         percent of juveniile killed by drug; could be 0 or avg micro/macro
#    mda_sterile : float
#         percent of adult worms perm sterilized
#    mda_clear : int
#         clearance time
#
#    Returns
#    -------
#    dfMF
#    dfAdult
#    dfJuv
#    dfHost
#    dfSel
#
#   '''
#    mda_start = mdalist[0]
#    mda_num = mdalist[1]
#    mda_freq = mdalist[2]
#    mda_coverage = mdalist[3]
#    mda_macro = mdalist[4]
#    mda_juvicide = mdalist[5]
#    mda_micro = mdalist[6]
#    mda_sterile = mdalist[7]
#    mda_clear = mdalist[8]
#
#    ##calculates time since mda as clear_count
#    if month < mda_start:
#            clear_count = 0
#    elif month >= mda_start and month < (mda_start + mda_freq * mda_num + mda_freq):
#        if month%mda_freq == 0:
#             clear_count = 1
#        elif (month - mda_clear)%mda_freq == 0:
#             clear_count = mda_clear + 1
#        else:
#             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
#                  clear_count = ((month-mda_freq)%mda_freq + 1)
#             else:
#                  clear_count = 0
#    elif month > (mda_start + mda_freq * mda_num):
#        clear_count = 0
###################################################3
#    ##assign MDA to hosts
#    if clear_count == 1:
#        ##add mda to dfHost per coverage
#        dfHost = hostmda_fx(village, dfHost, mda_coverage)
#         ##apply MDA effect to host populations of worms
#        for index, row in dfHost[dfHost.MDA == 1].iterrows():
#            #kill MF
#            try:
#                mdarandmf = np.random.random(len(dfMF.meta.hostidx == row.hostidx))
#                mdakillmf = np.where(mdarandmf < mda_micro ** dfMF.meta[dfMF.meta.hostidx == row.hostidx]["selS"])
#                dfMF.drop_worms(mdakillmf[0])
#            except ValueError:
#                pass
#             #kill Juv
#            try:
#                mdarandjuv = np.random.random(len(dfJuv.meta.hostidx == row.hostidx))
#                mdakilljuv = np.where(mdarandjuv < mda_juvicide ** dfJuv.meta[dfJuv.meta.hostidx == row.hostidx]["selS"])
#                dfJuv.drop_worms(mdakilljuv[0])
#            except ValueError:
#                pass
#            #kill Adult
#            try:
#                mdarandad = np.random.random(len(dfAdult.meta.hostidx == row.hostidx))
#                mdakillad = np.where(mdarandad < mda_macro ** dfAdult.meta[dfAdult.meta.hostidx == row.hostidx]["selS"])
#                dfAdult.drop_worms(mdakillad[0])
#            except ValueError:
#                pass
#########################################
#    ##normal surival funxtions
#    #adult worms and hosts are only evaluated per year
#    if month%12 == 0:
#        #Adult survival is based on weibull cdf
#        kill_adult_rand = np.random.random(dfAdult.meta.shape[0])
#        try:
#            kill_adultfxage = weibull_min.cdf(dfAdult.mata.age, shapeAdult,loc=0,scale=scaleAdult)
#        except TypeError:
#            kill_adultfxage = weibull_min.cdf(0, shapeAdult,loc=0,scale=scaleAdult)
#
#        kill_selfx = kill_adultfxage ** (1-abs(1-dfAdult.meta.selS))
#
#        dieAdult = np.where(kill_adult_rand < kill_selfx)
#        dfAdult.drop_worms(dieAdult[0])
#        dfAdult.meta.age = dfAdult.meta.age + 1
#        dfHost = dfHost[dfHost.age < dfHost.agedeath]
#        #remove all worms with dead host.hostidx from all dataframes
#        dfAdult.drop_worms(dfAdult.meta[~dfAdult.meta.hostidx.isin(dfHost.hostidx)].index.values)
#        dfJuv.drop_worms(dfJuv.meta[~dfJuv.meta.hostidx.isin(dfHost.hostidx)].index.values)
#        dfMF.drop_worms(dfMF.meta[~dfMF.meta.hostidx.isin(dfHost.hostidx)].index.values)
#        #add 1 year to all ages of hosts
#        dfHost.age = dfHost.age + 1
#        if hostmigrate != 0:
#            dfHost = hostmigration_fx(village, dfHost, hostmigrate)
#
#    ##Juv is exponential 0.866; surv_Juv
#    #dont include age 0 which just moved from transmission fx
#    dfJuv.age += 1
#    kill_juvrand = np.random.random(dfJuv.meta.shape[0])
#    kill_juvarray = np.repeat((1-surv_Juv), dfJuv.meta.shape[0])
#    kill_selfx = kill_juvarray ** (1-abs(1-dfJuv.meta.selS))
#    dieJuv = np.where(kill_juvrand < kill_selfx)
#    dfJuv.drop_worms(dieJuv[0])
#
#    ##MF is weibull cdf
#    kill_mfrand = np.random.random(dfMF.meta.shape[0])
#    try:
#        kill_mffxage = weibull_min.cdf(dfMF.meta.age,shapeMF,loc=0,scale=scaleMF)
#    except TypeError:
#        kill_mffxage = weibull_min.cdf(0,shapeMF,loc=0,scale=scaleMF)
#    kill_selfx = kill_mffxage ** (1-abs(1-dfMF.meta.selS))
#    dieMF = np.where(kill_mfrand < kill_selfx)
#    dfMF.drop_worms(dieMF[0])
#    dfMF.meta.age = dfMF.meta.age + 1 #2 - 12
#    dfMF.drop_worms(dfMF.meta.ix[dfMF.meta.age > 12].index.values)
#############################################
#    ##move Juv age 13 to adult age 1
#    juv_rows = dfJuv.meta[dfJuv.meta.age > 12].index.values
#    if any(juv_rows):
#        #reset age to adult
#        dfJuv.meta.ix[juv_rows,"age"] = 1
#        #increase R0net for next gen
#        dfJuv.meta.ix[juv_rows, "R0net"] += 1
#        dfAdult.add_worms(dfJuv, juv_rows)
#        dfJuv.drop_worms(juv_rows)
#
#    #fecundity calls mutation/recombination
#    dfAdult_mf, dfAdult = fecunditymda_sel2_fx(fecund, dfAdult, locus, mutation_rate,
#                                         recombination_rate, basepairs, selection,
#                                         densitydep_fec, mda_sterile, clear_count, mda_clear, dfHost)
#    dfAdult_mf.meta.sex = [random.choice("MF") for i in range(len(dfAdult_mf.meta))]
#    dfAdult_mf.meta.age = 1
#    dfMF.add_worms(dfAdult_mf, dfAdult_mf.meta.index.values)
#    return(dfHost, dfAdult, dfJuv, dfMF)