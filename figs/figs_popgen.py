#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import pandas as pd
from libsequence.polytable import simData
from libsequence.summstats import polySIM
from libsequence.windows import Windows
from libsequence.parallel import scheduler_init
init = scheduler_init(10)
from libsequence.fst import fst
#plotting
from .plottin.py import plot_allele_trace #for sel_trace_allele
from .plottin.py import plot_allele_frequency #for SFS and jSFS
from .plotting.py import plot_hapnetwork #for mitochondrial networks on MF
from .plotting.py import plot_pairwise

def site_freqspec_fx(dfworm):
    '''calculates the site frequency spectrum. If >1 pop also does jSFS

    Parameters
    ---------

    Returns
    -------
    '''


    plot_allele_frequency(SFS1)
    return(None)

def sel_trace_fx(dfworm):
    '''calculates the change in allele frequency through time for mutations with
    fitness affects

    Parameters
    ---------

    Returns
    -------
    '''

    plot_allele_trace(allele_trace)
    return(None)

def haplotype_net_fx(dfworm):
    '''constructs haplotype network table for plotting, vcf2hap

    Parameters
    ---------

    Returns
    -------
    '''
    mthaps = dfworm.h1['0']

    plot_hapnetwork(hapnetobj)
    return(None)

def figs2vcf_fx(dfworm):
    '''Takes worm obj and outputs vcf file

    Parameters
    ---------

    Returns
    -------
    '''


    return(None)

def figs2scikit_fx(dfworm):
    '''takes worm object and output format to load into scikit-allel

    Parameters
    ---------

    Returns
    -------
    '''


    return(None)

def dxy_fx(dfworm, mon):
    '''calculates pairwise sequence divergence between 2 pops

    Parameters
    ---------

    Returns
    -------
    '''
    dxy = ''
    write2file.csv("dxymatrix_{}".format(mon))
    plot_pairwise(dxy)
    return(None)

def thetaW_thetaB_fx(sdmf, mon):
    '''calculates the ratio of theta within hosts, to that of between hosts

    Parameters
    ---------

    Returns
    -------
    '''
    thetaWB = ''
    write2file.csv("TWTBmatrix_{}".format(mon))
    plot_pairwise(thetaWB)
    return(None)

def fst_pair_fx(sdmf, mon):
    '''calculates the pairwise FST between hosts

    Parameters
    ---------

    Returns
    -------
    '''
    #fst
    sdmf.size()
    f = fst(sdmf, [mfv1.shape[0], mfv2.shape[0]])
    f.hsm()
    write2file.csv("fstmatrix_{}".format(mon))
    plot_pairwise(fst)
    return(None)

def hostpopgen_fx(dfworm, outstats, mon):
    '''calculates popgen statistics for each host population

    Parameters
    ---------

    Returns
    -------
    '''
    sample_size = outstats[1]
    num_windows = outstats[3]
    basepairs = outstats[0]
    hostidx = dfworm.meta[(dfworm.meta.stage == "M")].hostidx.index.values
    for locus in range(1, len(basepairs)):
        win_len = float(outstats[2]) / basepairs[locus]
        win_step = 1.0/num_windows
        mfgeno = dfworm.meta[(dfworm.meta.stage == "M")].sample(sample_size).index.values
        #haplotype 1
        sdmf1 = simData()
        sdmf1.assign_sep(dfworm.pos[str(locus)] / basepairs[locus], dfworm.h1[str(locus)][mfgeno])
        psmf1 = polySIM(sdmf1,sizemf)
        #haplotype 2
        sdmf2 = simData()
        sdmf2.assign_sep(dfworm.pos[str(locus)] / basepairs[locus], dfworm.h2[str(locus)][mfgeno])
        psmf2 = polySIM(sdmf2,sizemf)
        #stats
        theta = psmf1.thetaw() + psmf2.thetaw()
        tajimasd =  psmf1.tajimasd() + psmf2.tajimasd()
        thetapi = psmf1.thetapi() + psmf2.thetapi()
        fulid = psmf1.fulid() + psmf2.fulid()
        fulidstar = psmf1.fulidstar() + psmf2.fulidstar()
        fulif = psmf1.fulif() + psmf2.fulif()
        fulifstar = psmf1.fulifstar() + psmf2.fulifstar()
        hprime = psmf1.hprime() + psmf2.hprime()
        numexternalmutations = psmf1.numexternalmutations() + psmf2.numexternalmutations()
        nummutations = psmf1.nummutations() + psmf2.nummutations()
        numpoly = psmf1.numpoly() + psmf2.numpoly()
        numsingletons = psmf1.numsingletons() + psmf2.numsingletons()
#        wallsb = psmf1.wallsb() + psmf2.wallsb()
#        wallsbprime = psmf1.wallsbprime() + psmf2.wallsbprime()
#        wallsq = psmf1.wallsq() + psmf2.wallsq()
#        lhaf = psmf1.lhaf() + psmf2.lhaf()
#        garudStats = psmf1.garudStats() + psmf2.garudStats()

    popgenTable = pd.DataFrame({"month" : [mon] * len(hostidx),
                                "hostidx" : hostidx,
                                "theta" : theta,
                                "tajD" : tajimasd,
                                "thetapi" : thetapi,
                                "fulid" : fulid,
                                "fulidstar" : fulidstar,
                                "fulif" : fulif,
                                "fulifstar" : fulifstar,
                                "hprime" : hprime,
                                "numexternalmutations" : numexternalmutations,
                                "nummutations" : nummutations,
                                "numpoly" : numpoly,
                                "numsingletons" : numsingletons
                                })
    popgenTable.to_csv()
    return(thetaHost)

def villpopgen_fx(dfworm, outstats, vill, mon):
    '''calculates popgen statistic for each village
    [basepairs, sample_size, window_length, num_windows, wb2vcf, figs2scikit]
    Parameters
    ---------
    dfworm : df
    outstats : list
    vill : int

    Returns
    -------

    '''
    sample_size = outstats[1]
    num_windows = outstats[3]
    basepairs = outstats[0]
    for locus in range(1, len(basepairs)):
        win_len = float(outstats[2]) / basepairs[locus]
        win_step = 1.0/num_windows

        adgeno = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "A")].sample(sample_size).index.values
        jvgeno = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "J")].sample(sample_size).index.values
        mfgeno = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "M")].sample(sample_size).index.values


        #''.join(str(n) for n in x)
        #[''.join(str(n) for n in y) for y in x]

        #haplotype 1
        sdad1 = simData()
        sdad1.assign_sep(dfworm.pos[str(locus)] / basepairs[locus], dfworm.h1[str(locus)][adgeno])
        sdjv1 = simData()
        sdjv1.assign_sep(dfworm.pos[str(locus)] / basepairs[locus], dfworm.h1[str(locus)][jvgeno])
        sdmf1 = simData()
        sdmf1.assign_sep(dfworm.pos[str(locus)] / basepairs[locus], dfworm.h1[str(locus)][mfgeno])

        psad1 = polySIM(sdad1,sizead)
        psjv1 = polySIM(sdjv1,sizejv)
        psmf1 = polySIM(sdmf1,sizemf)

        #haplotype 2
        sdad2 = simData()
        sdad2.assign_sep(dfworm.pos[str(locus)] / basepairs[locus], dfworm.h2[str(locus)][adgeno])
        sdjv2 = simData()
        sdjv2.assign_sep(dfworm.pos[str(locus)] / basepairs[locus], dfworm.h2[str(locus)][jvgeno])
        sdmf2 = simData()
        sdmf2.assign_sep(dfworm.pos[str(locus)] / basepairs[locus], dfworm.h2[str(locus)][mfgeno])

        psad2 = polySIM(sdad2,sizead)
        psjv2 = polySIM(sdjv2,sizejv)
        psmf2 = polySIM(sdmf2,sizemf)

        thetaA_t = psad1.thetaw() + psad2.thetaw()
        thetaJ_t = psjv1.thetaw() + psjv2.thetaw()
        thetaM_t = psmf1.thetaw() + psmf2.thetaw()
        tajD_t = psmf1.tajimasd() + psmf2.tajimasd()

        fst_b_t = fst_pair_fx(sdmf1, sdmf2, mon)
        thetaWB_t = thetaW_thetaB_fx(sdmf1, sdmf2, mon)
        dxy_t = dxy_fx(dfworm, mon)

    return(thetaA_t, thetaJ_t, thetaM_t, fst_b_t, thetaWB_t, dxy_t, tajD_t)
