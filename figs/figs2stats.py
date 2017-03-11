#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 19:29:56 2017

@author: stsmall
"""
import pickle
import numpy as np
import pandas as pd

#popgen
from .figs_popgen import villpopgen_fx #time series, only on MF
from .figs_popgen import figs2vcf_fx #only on final, all stages
from .figs_popgen import figs2scikit_fx
#demography
from .figs_demo import demo_stats_fx
from .figs_demo import prevTrans_fx
from .figs_demo import R0net_fx
from .figs_demo import host_stats_fx
from .figs_demo import demo_hoststats_fx

def output_tables_fx(logTime, sim_time, outstats):
    '''builds summary table from FiGS simulations

    Parameters
    ---------

    Returns
    -------
    '''

    #prev, transmission
    with open('L3transdict.pkl','rb') as input:
        L3transdict = pickle.load(input)
    avgTrans, avgPrev, varTrans, varPrev, villages = prevTrans_fx(L3transdict, logTime)

    #R0net
    with open('R0netlist.pkl','rb') as input:
        R0netlist = pickle.load(input)
    R0, ravg, rvar = R0net_fx(R0netlist)

    #demo table params
    month = np.repeat(range(logTime, sim_time, logTime), villages)
    village = range(villages) * len(month)
    infhost = []
    adult = []
    juv = []
    mf = []
    vadult = []
    vjuv = []
    vmf = []
    thetaA = []
    thetaJ = []
    thetaM = []
    fst_b = []
    da =[]
    dxy = []
    tajD =[]
    pi = []

    ####load data from incremental output
    for mon in range(logTime, sim_time, logTime):
        print("figs2stats")
        #worm files
        with open('dfworm_{}.pkl'.format(mon), 'rb') as worm:
            dfworm = pickle.load(worm)

        for vill in range(villages):
            #demo stats
            adult_t, vadult_t, juv_t, vjuv_t, mf_t, vmf_t = demo_stats_fx(dfworm, vill)

            adult.append(adult_t)
            juv.append(juv_t)
            mf.append(mf_t)
            vadult.append(vadult_t)
            vjuv.append(vjuv_t)
            vmf.append(vmf_t)

            #popgen stats
            thetaA_t, thetaJ_t, thetaM_t, fst_b_t, da_t, dxy_t, tajD_t, pi_t= villpopgen_fx(dfworm, outstats, vill, mon)

            thetaA.append(thetaA_t)
            thetaJ.append(thetaJ_t)
            thetaM.append(thetaM_t)
            fst_b.append(fst_b_t)
            da.append(da_t)
            dxy.append(dxy_t)
            tajD.append(tajD_t)
            pi.append(pi_t)

        ##host stats
        with open('dfHost_{}.pkl'.format(mon), 'rb') as host:
            dfHost = pickle.load(host)
        infhost_t = host_stats_fx(dfHost, villages)
        infhost.append(infhost_t)
        demo_hoststats_fx(dfworm, dfHost, mon)

    #demotable
    summaryTable = pd.DataFrame({"month" : month,
                                "village" : village,
                                "inf_host" : infhost,
                                "avg_prev" : avgPrev,
                                "var_prev" : varPrev,
                                "avg_trans" : avgTrans,
                                "var_trans" : varTrans,
                                "R0" : R0,
                                "avg_repo" : ravg,
                                "var_repo" : rvar,
                                "avg_adult" : adult,
                                "avg_juv" : juv,
                                "avg_mf" : mf,
                                "var_adult" : vadult,
                                "var_juv" : vjuv,
                                "var_mf" : vmf,
                                "thetaA" : thetaA,
                                "thetaJ" : thetaJ,
                                "thetaM" : thetaM,
                                "fst_b" : fst_b,
                                "da" : da,
                                "dxy" : dxy,
                                "tajD" : tajD,
                                "pi" : pi
                                })
    summaryTable = summaryTable.loc[:, ['month', 'village', 'inf_host', 'avg_prev', 'var_prev', 'avg_trans',
                                  'var_trans', 'R0', 'avg_repo', 'var_repo', 'avg_adult', 'avg_juv', 'avg_mf',
                                  'var_adult', 'var_juv', 'var_mf','thetaA', 'thetaJ', 'thetaM', 'fst_b', 'da', 'dxy', 'tajD', 'pi']]
    summaryTable.to_csv(summaryTable)


    #####final dataframe from figs
    if outstats[4] or outstats[5]:
        with open('dfworm_final.pkl', 'rb') as worm:
            dfworm = pickle.load(worm)
        if outstats[4]:
            figs2vcf_fx(dfworm)
        else: pass
        if outstats[5]:
            figs2scikit_fx(dfworm)
        else: pass
    else: pass

    return(None)
