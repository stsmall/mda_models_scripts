#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 19:29:56 2017

@author: stsmall
"""
import pickle
import numpy as np
import pandas as pd

def figs2pylibseq_fx(dfworm, outstats):
    '''

    '''
    #outstats = [sample_size, window_length, num_windows, wb2vcf, wbmfvcf, wbadultvcf,
    #            wbjuvvcf, wbfracvcf, figs2scikit]
    with open('dfworm_60.pkl','rb') as input:
        dfworm = pickle.load(input)

    sample_size = outstats[0]
    window_length = outstats[1]
    num_windows = outstats[2]

    return None

def figs2scikit_fx(dfworm, outstats):
    '''

    '''
    #outstats = [sample_size, window_length, num_windows, wb2vcf, wbmfvcf, wbadultvcf,
    #            wbjuvvcf, wbfracvcf, figs2scikit]
    with open('dfworm_60.pkl','rb') as input:
        dfworm = pickle.load(input)


    return None

def R0net_fx(R0netlist):
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
    with open('R0netlist.pkl','rb') as input:
        R0netlist = pickle.load(input)

    R0 = [float(j) / R0netlist[i-1] for i, j in enumerate(R0netlist)][1:]

    return(R0)

def prevTrans_fx(R0netlist):
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
    with open('L3transdict.pkl','rb') as input:
        L3transdict = pickle.load(input)

    trans = [L3transdict[str(vill)][0] for vill in range(len(L3transdict.keys()))]
    prev = [L3transdict[str(vill)][1] for vill in range(len(L3transdict.keys()))]
    return(trans, prev)