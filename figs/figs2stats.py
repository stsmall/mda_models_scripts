#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 19:29:56 2017

@author: stsmall
"""
import pickle
import numpy as np
import pandas as pd
from libsequence.polytable import simData
from libsequence.summstats import polySIM
from libsequence.windows import simDataWindows
from libsequence.fst import fst
from .plottin.py import plot_coordinates_host
from .plottin.py import sel_allele_trace
from .plottin.py import plot_allele_frequency

def figs2pylib_fx(dfworm, outstats):
    '''
    population genentic stats
    '''
    #outstats = [sample_size, window_length, num_windows, wb2vcf, wbmfvcf, wbadultvcf,
    #            wbjuvvcf, wbfracvcf, figs2scikit]

    sample_size = outstats[0]
    mfv1 = dfworm.meta[('village == 1') & ('stage == "M"')].sample(frac = sample_size).index.values
    mfv2 = dfworm.meta[('village == 2') & ('stage == "M"')].sample(frac = sample_size).index.values
#    jv1 =
#    jv2 =
#    av1 =
#    av2 =

    mfgeno = np.append(mfv1,mfv2)
    basepairs = [13000]
    #pylibseq summary of locus
    sdmf = simData()
    sdmf.assign_sep(dfworm.pos['0'] / basepairs[0], dfworm.h1['0'][mfgeno])
    ps = polySIM(sdmf)
    theta = ps.thetaw()
    tajimaD = ps.tajimasd()

    #windows
    window_length = outstats[1]
    num_windows = outstats[2]
    size = window_length / basepairs[0]
    step = 1.0 / num_windows
    w = simDataWindows(sdmf, window_size=size, step_len=step, starting_pos=0., ending_pos=1.0)
    for i in range(len(w)):
        wi = w[i]
        pswi = polySIM(wi)
        print(pswi.thetaw())

    #fst
    sdmf.size()
    f = fst(sdmf, [mfv1.shape[0], mfv2.shape[0]])
    f.hsm()

    return(theta, tajimaD)

def figs2scikit_fx(dfworm, outstats):
    '''

    '''
    #outstats = [sample_size, window_length, num_windows, wb2vcf, wbmfvcf, wbadultvcf,
    #            wbjuvvcf, wbfracvcf, figs2scikit]


    return None

def popgen_stats_fx(dfworm, outstats):
    '''
    '''
    theta, tajimaD = figs2pylib_fx(dfworm, outstats)
    #hapdiv
    #hapdiv = np.unique(dfworm.h1['0'],return_counts=True) / dfworm.h1['0'].shape[0]
    hapdiv = ''
    #sel_af
    sel_af = ''

    return(theta, tajimaD, hapdiv, sel_af)
def prevTrans_fx(L3transdict, logTime):
    '''Calculates the reproductive number, R0, by counting the uniqueness
    of R0net per village and taking the mean counts

    Parameters
    ----------
    dfJuv : df
        dataframe of juvenilles age 13

    Returns
    --------
    R0 : float, list
        reproductive rate of each village
    '''
    #logTime = 12
    with open('L3transdict.pkl','rb') as input:
        L3transdict = pickle.load(input)

    #trans_prev
    vill = len(L3transdict['prev'][0])
    time = len(L3transdict['prev'])
    trans = [e for l in L3transdict['trans'] for e in l]
    prev = np.round([e for l in L3transdict['prev'] for e in l],2)
    month = np.repeat(range(1, time + 1), vill)
    village = range(vill) * time
    transTable = pd.DataFrame({"month" : month,
                                "village" : village,
                                "trans" : trans,
                                "prev" : prev})
    transTable = transTable.loc[:, ['month','village', 'trans', 'prev']]
    transTable.to_csv(transTable)

    meantran = [transTable.groupby("village").rolling(logTime).mean().dropna()['trans'][v::logTime][v] for v in range(vill)]
    vartran = [transTable.groupby("village").rolling(logTime).mean().dropna()['trans'][v::logTime][v] for v in range(vill)]
    meanprev = [transTable.groupby("village").rolling(logTime).mean().dropna()['prev'][v::logTime][v] for v in range(vill)]
    varprev = [transTable.groupby("village").rolling(logTime).var().dropna()['prev'][v::logTime][v] for v in range(vill)]
    avgTrans = [val for pair in zip(*meantran) for val in pair]
    varTrans = [val for pair in zip(*vartran) for val in pair]
    avgPrev = [val for pair in zip(*meanprev) for val in pair]
    varPrev = [val for pair in zip(*varprev) for val in pair]
    return(avgTrans, avgPrev, varTrans, varPrev, vill)

def demostats_fx(logTime, sim_time, outstats):
    '''Calculates the reproductive number, R0, by counting the uniqueness
    of R0net per village and taking the mean counts

    Parameters
    ----------
    dfJuv : df
        dataframe of juvenilles age 13

    Returns
    --------
    R0 : float, list
        reproductive rate of each village
    '''
    #prev, transmission
    with open('L3transdict.pkl','rb') as input:
        L3transdict = pickle.load(input)
    #returns avg of 12 months
    avgTrans, avgPrev, varTrans, varPrev, vill= prevTrans_fx(L3transdict, logTime)

    #R0net
    with open('R0netlist.pkl','rb') as input:
        R0netlist = pickle.load(input)
    #R0 by village
    R0 = []
    R0_t = zip(*R0netlist['R0'])
    for v in range(len(R0_t)):
        R0.append([float(j) / R0_t[v][i-1] for i,j in enumerate(R0_t[v])])

    #reproductive mean and var
    ravg = [j for time in R0netlist['repoavg'] for j in time]
    rvar = [j for time in R0netlist['repovar'] for j in time]

    #demo table params
    month = np.repeat(range(logTime, sim_time, logTime), vill)
    village = range(vill) * len(month)
    infhost = []
    adult = []
    juv = []
    mf = []
    vadult = []
    vjuv = []
    vmf = []
    theta = []
    tajD = []
    hapdiv = []
    sel_allelefreq = []

    #load data
    for mon in range(logTime, sim_time, logTime):
        #worm stats
        with open('dfworm_{}.pkl'.format(mon), 'rb') as worm:
            dfworm = pickle.load(worm)
        for v in range(vill):
            adult.append(dfworm.meta.groupby(["stage","village","hostidx"]).size()['A'][v].mean())
            vadult.append(dfworm.meta.groupby(["stage","village","hostidx"]).size()['A'][v].var())
            juv.append(dfworm.meta.groupby(["stage","village","hostidx"]).size()['J'][v].mean())
            vjuv.append(dfworm.meta.groupby(["stage","village","hostidx"]).size()['J'][v].var())
            mf.append(dfworm.meta.groupby(["stage","village","hostidx"]).size()['M'][v].mean())
            vmf.append(dfworm.meta.groupby(["stage","village","hostidx"]).size()['M'][v].var())

        #popgen stats
        theta_v, tajD_v, hapdiv_v, sel_af = popgen_stats_fx(dfworm, outstats)
        theta.append(theta_v)
        tajD.append(tajD_v)
        hapdiv.append(hapdiv_v)
        sel_allelefreq.append(sel_af)

        #host stats
        with open('dfHost_{}.pkl'.format(mon), 'rb') as host:
            dfHost = pickle.load(host)
        infhost.append([dfHost.groupy("village").size()[i] for i in range(len())])

    ##final dataframe for figs
    with open('dfworm_final.pkl', 'rb') as worm:
        dfworm = pickle.load(worm)
    plot_allele_frequency(dfworm)
    figs2scikit_fx(dfworm, outstats)
    sel_allele_trace(sel_allelefreq)

    with open('dfHost_final.pkl','rb') as host:
        dfHost = pickle.load(host)
    plot_coordinates_host(dfHost)

    #demotable
    demoTable = pd.DataFrame({"month" : month,
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
                                "theta" : theta,
                                "tajD" : tajD,
                                "hapdiv" : hapdiv
                                })
    demoTable = demoTable.loc[:, ['month', 'village', 'inf_host', 'avg_prev', 'var_prev', 'avg_trans',
                                  'var_trans', 'R0', 'avg_repo', 'var_repo', 'avg_adult', 'avg_juv', 'avg_mf',
                                  'var_adult', 'var_juv', 'var_mf','theta', 'tajD', 'hapdiv']]
    demoTable.to_csv(demoTable)

    return(None)

#if __name__ == '__main__':
#    demostats_fx(logTime, sim_time, outstats)
