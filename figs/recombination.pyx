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
from libcpp.list cimport list as cpplist
#from cython.parallel import parallel, prange
#from libc.stdlib cimport abort, malloc, free

DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t



@cython.boundscheck(False)
@cython.wraparound(False)
cdef long[: ,::1] mate_worms(
        long[:] mate_array, 
        long[:] fec,
        np.ndarray[np.uint64_t, ndim=1] pos,
        int basepairs,
        float recomb_rate,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] fh1, 
        np.ndarray[DTYPE_t, ndim=2, mode='c'] fh2,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] mh1,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] mh2):
    """ Mates and recombines at a given loci
    """
    # :TODO need to check max integer
    cdef np.intp_t i, outsize 
    outssize = np.sum(fec) 
    cdef np.ndarray mnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    cdef np.ndarray fnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    # Haplotype chooser
    cdef np.ndarray hapc = np.random.randit(0, 2,
            outsize)
    h1 = np.zeros((outssize, fh1.shape[1]), dtype=np.uint8)
    h2 = np.zeros((outssize, fh1.shape[1]), dtype=np.uint8)
    for i in range(outsize):
        if mnum_recomb[i] == 0:
            h1[i, :] = fh1[i, :]
        elif fnum_recomb[i] == 0:
            pass
        else: 
            pass
    bleh = np.empty( (20, 20), dtype=np.int)
    return(bleh)


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
    dfAdult_mf : figs.worm.Worms object
        Worms containing new larval parasites
    dfAdult : figs.worm.Worms object
        dataframe containing reproducing adults
    recombination_rate : float, list
        recombination rate for each locus
    basepairs : int, list
        length of each locus in basepairs

    Returns
    -------
    dfAdult_mf : pd df

    """
    hosts = dfAdult.meta.hostidx.unique()
    cdef str host
    # How to type this?
    #cdef bool[:] ahost, females, males
    cdef Py_ssize_t loc
    cdef cpplist[int] boundry_list
    cdef long[:] fec 
    #cdef int[:] mate_array = np.empty(np.sum(females))
    for loc in range(locus):
        if recombination_rate[loc] != 0:
            boundry_list.push_back(dfAdult.h1[str(loc)].shape[1])
        else: pass
    for host in hosts:
        #chost = dfAdult.meta.index[dfAdult.meta.hostidx == host].values
        ahost = dfAdult.meta.hostidx == host
        females = np.logical_and(ahost, dfAdult.meta.sex == 'F').values
        males = np.logical_and(ahost, dfAdult.meta.sex == 'M').values
        if np.sum(males) == 0 or np.sum(females) == 0:
            continue
        else:
            fec = dfAdult.meta.fec[females].values
            mate_array = np.random.randint(0, np.sum(males), np.sum(females),
                    dtype=np.int)
            for loc in range(locus):
                if recombination_rate[loc] == 0:
                    pass
                else:
                    print(mate_array)
                    mate_worms(mate_array, 
                            fec, 
                            dfAdult.pos[str(loc)],
                            basepairs[loc],
                            recombination_rate[loc],
                            dfAdult.h1[str(loc)][females, :],
                            dfAdult.h2[str(loc)][females, :],
                            dfAdult.h1[str(loc)][males, :],
                            dfAdult.h2[str(loc)][males, :])
    '''
    cdef int N
    cdef int i
    lid = "locus_{0!s}"
    dout = []
    loci = [Locus(lid.format(i), recombination_rate = j, basepairs = k)
            for i, j, k in zip(range(locus), recombination_rate, basepairs) if
            j != 0]
    hosts = dfAdult.hostidx.unique()
    N = hosts.shape[0]
    for i in range(N):
        dout.append(_temp(dfAdult[dfAdult.hostidx == hosts[i]], loci))
    dfAdult_mf = pd.concat(dout)
    dfAdult_mf.reset_index(drop=True, inplace=True)
    '''
    return(dfAdult)
