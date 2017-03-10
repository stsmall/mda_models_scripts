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
from libsequence.summstats import garudStats
from libsequence.summstats import ld
from libsequence.summstats import lhaf
from libsequence.windows import Windows
from libsequence.parallel import scheduler_init
init = scheduler_init(10)
from libsequence.fst import fst
#plotting
from .plottin.py import plot_allele_trace #for sel_trace_allele
from .plottin.py import plot_allele_frequency #for SFS and jSFS
from .plotting.py import plot_hapnetwork #for mitochondrial networks on MF
from .plotting.py import plot_pairwise

def site_freqspec_fx(dfworm, mon, mf, locus):
    '''calculates the site frequency spectrum. If >1 pop also does jSFS

    Parameters
    ---------

    Returns
    -------
    '''
    #calculate SFS
    mfsfs1 = dfworm.h1[locus][mf]
    mfsfs2 = dfworm.h2[locus][mf]
    freqsum = sum(dfworm.h1[locus][mfsfs1]) + sum(dfworm.h2[locus][mfsfs2])
    sfs = np.unique(np.sort(freqsum), return_counts = True)[1]
    plot_allele_frequency(sfs)

    #pass to sel_trace
    sel_trace_fx(dfworm, freqsum, locus, mf)
    return(None)

def sel_trace_fx(dfworm, freqsum, locus, mf):
    '''calculates the change in allele frequency through time for mutations with
    fitness affects

    Parameters
    ---------

    Returns
    -------
    '''
    if bool(dfworm.sel):
        #records trace of selected alleles
        ftrace = np.where(dfworm.sel[locus + 'F'] > 0)[0]
        f_freq = freqsum[ftrace] / (2.0 * mf.shape[0])
        strace = np.where(dfworm.sel[locus + 'S'] > 0)[0]
        s_freq = freqsum[strace] / (2.0 * mf.shape[0])
    plot_allele_trace(f_freq, s_freq)
    return(None)

def haplotype_net_fx(dfworm):
    '''constructs haplotype network table for plotting, vcf2hap

    Parameters
    ---------

    Returns
    -------
    '''
    #mthaps = dfworm.h1['0']
    #vcf2networks, fitchi, pegas/adgenet in R ???

    #plot_hapnetwork(hapnetobj)
    return(None)

def figs2vcf_fx(dfworm):
    '''Takes worm obj and outputs vcf file

    Parameters
    ---------

    Returns
    -------
    '''


    return(None)

def figs2scikit_fx(dfworm, sample_size, vill):
    '''takes worm object and output format to load into scikit-allel

    Parameters
    ---------

    Returns
    -------
    '''
    import allel
    sample_size = 30
    vill = 1
    #groupby: village == vill, stage == "M", hostidx
    mf = dfworm.meta[(dfworm.meta["stage"] == "M") & (dfworm.meta["village"] == vill)].index.values
    mf_subsample = mf[np.random.choice(mf.shape[0], sample_size, replace = False)]
    haps1 = dfworm.h1['1'][mf_subsample]
    haps2 = dfworm.h2['1'][mf_subsample]
    #import into allel
    positions = dfworm.pos['1']
    haparray_ind = np.dstack((haps1, haps2))
    haparray_pop = np.concatenate(haparray_ind, axis = 1)
    h = allel.HaplotypeArray(haparray_pop)
    #g = h.to_genotypes(ploidy=2)
    ac = h.count_alleles()
#    mf_mask = dfworm.meta.ix[mf].groupby("hostidx").size() < sample_size
#    if any(mf_mask):
#        mf = dfworm.meta.ix[mf][~dfworm.meta["hostidx"].isin(mf_mask[mf_mask].index.values)].index.values
#    mf_pairs = dfworm.meta.ix[mf].groupby("hostidx").apply(lambda y: y.sample(sample_size).index.values)
#    hostidx = [h for h in mf_pairs.keys()]

    return(None)

def pairwise_div_fx(dfworm, mon, vill, basepairs, sample_size):
    '''calculates the pairwise FST between hosts

    Parameters
    ---------

    Returns
    -------
    fst_slat : float
        slatkin 1991 calculation of Fst, mean over all host-pops in village
    '''
    fst_t = []
    dxy = []
    da = []
    theta = []
    tajimasd = []
    thetapi = []
    fulid = []
    fulidstar = []
    fulif = []
    fulifstar = []
    hprime = []
    numexternalmutations = []
    nummutations = []
    numpoly = []
    numsingletons = []
    pi = []
#####
    #groupby: village == vill, stage == "M", hostidx
    mf = dfworm.meta[(dfworm.meta["stage"] == "M") & (dfworm.meta["village"] == vill)].index.values
    mf_mask = dfworm.meta.ix[mf].groupby("hostidx").size() < sample_size
    if any(mf_mask):
        mf = dfworm.meta.ix[mf][~dfworm.meta["hostidx"].isin(mf_mask[mf_mask].index.values)].index.values
    mf_pairs = dfworm.meta.ix[mf].groupby("hostidx").apply(lambda y: y.sample(sample_size).index.values)
    hostidx = [h for h in mf_pairs.keys()]
#####
    for loc in range(1, len(basepairs)):
        locus = str(loc)
        #print sfs and allele traces
        site_freqspec_fx(dfworm, mon, mf, locus)
        for host1 in hostidx:
            for host2 in hostidx:
                if host1 != host2:
                    pos = dfworm.pos[locus] / float(basepairs[loc])
                    #pylibseq
                    pop1 = np.vstack([dfworm.h1[locus][mf_pairs[host1]], dfworm.h2[locus][mf_pairs[host1]]])
                    pop2 = np.vstack([dfworm.h1[locus][mf_pairs[host2]], dfworm.h2[locus][mf_pairs[host2]]])
                    gtpop1 = [''.join(str(n) for n in y) for y in pop1]
                    gtpop2 = [''.join(str(n) for n in y) for y in pop2]

                    #host stats
                    sdpop1 = simData()
                    sdpop1.assign_sep(pos, gtpop1)
                    pspop1 = polySIM(sdpop1)
                    sdpop2 = simData()
                    sdpop2.assign_sep(pos, gtpop2)
                    pspop2 = polySIM(sdpop2)

                    theta.append(pspop1.thetaw())
                    theta.append(pspop2.thetaw())
                    tajimasd.append(pspop1.tajimasd())
                    tajimasd.append(pspop2.tajimasd())
                    thetapi.append(pspop1.thetapi())
                    thetapi.append(pspop2.thetapi())
                    fulid.append(pspop1.fulid())
                    fulid.append(pspop2.fulid())
                    fulidstar.append(pspop1.fulidstar())
                    fulidstar.append(pspop2.fulidstar())
                    fulif.append(pspop1.fulif())
                    fulif.append(pspop2.fulif())
                    fulifstar.append(pspop1.fulifstar())
                    fulifstar.append(pspop2.fulifstar())
                    hprime.append(pspop1.hprime())
                    hprime.append(pspop2.hprime())
                    numexternalmutations.append(pspop1.numexternalmutations())
                    numexternalmutations.append(pspop2.numexternalmutations())
                    nummutations.append(pspop1.nummutations())
                    nummutations.append(pspop2.nummutations())
                    numpoly.append(pspop1.numpoly())
                    numpoly.append(pspop2.numpoly())
                    numsingletons.append(pspop1.numsingletons())
                    numsingletons.append(pspop2.numsingletons())

#                    wallsb = pspop1.wallsb()
#                    wallsbprime = pspop1.wallsbprime()
#                    wallsq = pspop1.wallsq()
#                    garudStats_t = garudStats(sdpop1)
#                    lhaf_t = lhaf(sdpop1,10) #what is the double?
#                    ld_t = ld(sdpop1, haveOutgroup = False, mincount = .05, maxDist = 5000)

                    #fst
                    sdfst = simData()
                    geno_fst = gtpop1 + gtpop2
                    sdfst.assign_sep(pos, geno_fst)
                    size = [pop1.shape[0], pop2.shape[0]]
                    f1 = fst(sdfst, size)
                    fst_t.append(f1.slatkin())

                    #dxy, sample sizes must be equal
                    sizedxy = min(size)
                    dxy_t = sum([sum((i + j) == 1) for i in pop1[:sizedxy] for j in pop2[:sizedxy]]) / (sizedxy**2)
                    pi_p1 = sum([sum((x + y) == 1) for i, x in enumerate(pop1) for j, y in enumerate(pop1) if i != j ]) / (len(pop1) * (len(pop1) - 1 ) / 2.0)
                    pi_p2 = sum([sum((x + y) == 1) for i, x in enumerate(pop2) for j, y in enumerate(pop2) if i != j ]) / (len(pop2) * (len(pop2) - 1 ) / 2.0)
                    da.append(dxy_t - ((pi_p1 + pi_p2) / 2.0))
                    dxy.append(dxy_t)
                    pi.append(pi_p1)
                    pi.append(pi_p2)

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
    popgenTable.to_csv("popgenHostTable_v{}_m{}".format(vill, mon))

#    #writes full matrix for analysis
#    write2file.csv("fstmatrix_v{}_m{}".format(vill, mon))
#    write2file.csv("dxymatrix_v{}_m{}".format(vill, mon))
#    write2file.csv("damatrix_v{}_m{}".format(vill, mon))
#    #surface plot of pairwise
#    plot_pairwise(fst)
#    plot_pairwise(dxy)
#    plot_pairwise(da)
    return(fst, dxy, da, pi)

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
    #vary these for sample size and locus size info
    sample_size = outstats[1] #number of sequences to sample
    num_windows = outstats[3]
    basepairs = outstats[0]

    ##example of 1 window of 5kb from the start of a 13kb locus
    #window_length = outstats[2] #5000
    #win_size = window_length / basepairs[loc] #5000 / 13000.0
    #num_windows = 1
    #win_step = 1.0 / num_windows #1.0
    #w = Windows(sd, window_size = 0.38, step_len = 1, starting_pos = 0., ending_pos = 1.0)
    #w = Windows(sd, window_size = 0.38, step_len = 0.38, starting_pos = 0.38, ending_pos = 0.76)

    thetaA_t = []
    thetaJ_t = []
    thetaM_t = []
    tajD_t = []

    for loc in range(1, len(basepairs)):
        locus = str(loc)
        pos = dfworm.pos[locus] / float(basepairs[loc])
#        win_size = float(outstats[2]) / float(basepairs[loc])
#        win_step = 1.0/num_windows
        adpop = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "A")].sample(sample_size).index.values
        jvpop = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "J")].sample(sample_size).index.values
        mfpop = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "M")].sample(sample_size).index.values

        mfgeno1 = [''.join(str(n) for n in y) for y in dfworm.h1[locus][mfpop]]
        mfgeno2 = [''.join(str(n) for n in y) for y in dfworm.h2[locus][mfpop]]
        jvgeno1 = [''.join(str(n) for n in y) for y in dfworm.h1[locus][jvpop]]
        jvgeno2 = [''.join(str(n) for n in y) for y in dfworm.h2[locus][jvpop]]
        adgeno1 = [''.join(str(n) for n in y) for y in dfworm.h1[locus][adpop]]
        adgeno2 = [''.join(str(n) for n in y) for y in dfworm.h2[locus][adpop]]

        #haplotype 1
        sdad1 = simData()
        sdad1.assign_sep(pos, adgeno1)
        sdjv1 = simData()
        sdjv1.assign_sep(pos, jvgeno1)
        sdmf1 = simData()
        sdmf1.assign_sep(pos, mfgeno1)
        psad1 = polySIM(sdad1)
        psjv1 = polySIM(sdjv1)
        psmf1 = polySIM(sdmf1)

        #haplotype 2
        sdad2 = simData()
        sdad2.assign_sep(pos, adgeno2)
        sdjv2 = simData()
        sdjv2.assign_sep(pos, jvgeno2)
        sdmf2 = simData()
        sdmf2.assign_sep(pos, mfgeno2)
        psad2 = polySIM(sdad2)
        psjv2 = polySIM(sdjv2)
        psmf2 = polySIM(sdmf2)

        #stats
        thetaA_t.append((psad1.thetaw() + psad2.thetaw()) / 2.0)
        thetaJ_t.append((psjv1.thetaw() + psjv2.thetaw()) / 2.0)
        thetaM_t.append((psmf1.thetaw() + psmf2.thetaw()) / 2.0)
        tajD_t.append((psmf1.tajimasd() + psmf2.tajimasd()) / 2.0)

    #call paired functions
    fst, dxy, da, pi = pairwise_div_fx(dfworm, mon, vill, basepairs, sample_size) #returns mean fst_slatkin
    fst_slat_t = np.mean(fst)
    dxy_t = np.mean(dxy)
    da_t = np.mean(da)
    pi_t = np.mean(pi)

    return(np.mean(thetaA_t), np.mean(thetaJ_t), np.mean(thetaM_t), fst_slat_t,
           da_t, dxy_t, tajD_t, pi_t)