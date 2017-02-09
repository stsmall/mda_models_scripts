#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
FiGS Copyright (C) 2017 Scott T. Small
This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
This is free software, and you are welcome to redistribute it
under certain conditions; type `show c' for details.
'''
import numpy as np
cimport numpy as np
import random
import pandas as pd
import cython
from cython.parallel import parallel, prange
from libc.stdlib cimport abort, malloc, free


DTYPE = np.int
ctypedef np.int_t DTYPE_t

from .locus import Locus


def recombination_locus(np.ndarray[int, ndim=1] h1, 
        np.ndarray[int, ndim=1] h2, 
        int crossover_pos):
    """Calculates the recombination at a given locus
    """
    cdef int h1_ix = h1.shape[0]
    cdef int h2_ix = h2.shape[0] 
    cdef unsigned i = 0
    cdef int j
    for j in h1:
        if j > crossover_pos:
            h1_ix = i
            break
        else:pass
        i += 1
    i = 0
    for j in h2:
        if j > crossover_pos:
            h2_ix = i
            break
        else:pass
        i += 1
    h1_new = np.append(h1[0:h1_ix], h2[h2_ix:])
    h2_new = np.append(h2[0:h2_ix], h1[h1_ix:])
    return(h1_new, h2_new)


def _temp(df, loci):
    """ 
    Parameter
    ---------
    df : pandas dataframe
    loci : figs.Locus list

    Returns
    -------
    """
    cdef int num_recomb

    females = df.query('sex == "F" and fec > 0')
    males = df.ix[df.sex == "M", :]
    
    if males.shape[0] == 0:
        return(males)
    elif females.shape[0] == 0:
        return(females)
    outs = []
    for _, f_row in females.iterrows():
        nr = f_row.copy()
        male = males.sample(1).iloc[0, :]
        for mf in range(f_row.fec):
            for loc in loci:
                lid = loc.idx + '_h{0!s}'
                num_recomb = np.random.poisson(
                        loc.recombination_rate * loc.basepairs * 2)

                if num_recomb == 0:
                    nr[lid.format(1)] = f_row[lid.format(random.choice("12"))]
                    nr[lid.format(2)] = male[lid.format(random.choice("12"))]
                else:
                    sex_xing = random.choice("MF")
                    h1m = male[lid.format(1)]
                    h2m = male[lid.format(2)]
                    h1f = f_row[lid.format(1)]
                    h2f = f_row[lid.format(2)]
                    for _ in range(num_recomb):
                        crossover_pos = random.randint(0,
                                loc.basepairs)
                        if sex_xing == "M":
                            h1m, h2m = recombination_locus(h1m, h2m,
                                    crossover_pos)
                        elif sex_xing == "F":
                            h1f, h2f = recombination_locus(h1f, h2f,
                                    crossover_pos)
                    nr[lid.format(1)] = random.choice([h1f, h2f])
                    nr[lid.format(2)] = random.choice([h1m, h2m])
            outs.append(nr)
    outdf = pd.concat(outs, axis=1, ignore_index=True).T
    return(outdf)


@cython.boundscheck(False)
def recombination_fx(locus,
                     dfAdult,
                     list recombination_rate,
                     list basepairs):
    """calculate number of recombination events and rearranges haplotypes
    :TODO add for recombination map

    Parameters
    ---------
    locus: int
        number of loci
    dfAdult_mf : pandas dataframe
        dataframe containing new larval parasites
    dfAdult : pd df
        dataframe containing reproducing adults
    recombination_rate : float, list
        recombination rate for each locus
    basepairs : int, list
        length of each locus in basepairs

    Returns
    -------
    dfAdult_mf : pd df

    """
    cdef int N
    cdef int i
    lid = "locus_{0!s}"
    dout = []
    loci = [Locus(lid.format(i), recombination_rate = j, basepairs = k)
            for i, j, k in zip(range(locus), recombination_rate, basepairs) if 
            j != 0]
    hosts = dfAdult.hostidx.unique().value 
    N = hosts.shape[0]
    for i in range(N):
        dout.append(_temp(hosts[i], loci))
    dfAdult_mf = pd.concat(dout)
    return dfAdult_mf
