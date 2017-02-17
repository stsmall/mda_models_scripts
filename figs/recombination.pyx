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
from libc.stdlib cimport rand, RAND_MAX
#from cython.parallel import parallel, prange
#from libc.stdlib cimport abort, malloc, free

DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
cdef long[: ,::1] mate_worms(
        long[:] mate_array, 
        np.ndarray[np.int64_t, ndim=1] fec,
        np.ndarray[np.uint64_t, ndim=1] pos,
        int basepairs,
        float recomb_rate,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] fem, 
        np.ndarray[DTYPE_t, ndim=2, mode='c'] males):
    """ Mates and recombines at a given loci
    """
    # :TODO need to check max integer
    cdef np.intp_t i, outsize 
    outssize = np.sum(fec) 
    cdef int mnworms, fnworms
    cdef int mhapc, fhapc 
    cdef np.ndarray iix_ma = np.repeat(mate_array, 
            fec)
    cdef np.ndarray femindex = np.arange(fem.shape[0]/2)
    cdef np.ndarray iix_fem = np.repeat(femindex, fec)
    cdef np.ndarray mnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    cdef np.ndarray fnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    # Haplotype chooser
    mnworms = males.shape[0]/2
    fnworms = fem.shape[0]/2
    h1 = np.zeros((outssize, fem.shape[1]), dtype=np.uint8)
    h2 = np.zeros((outssize, fem.shape[1]), dtype=np.uint8)
    for i in range(outsize):
        if mnum_recomb[i] == 0:
            mhapc = int(rand()/RAND_MAX)
            print(str(mhapc))
            h1[i, :] = males[iix_ma[i] + mnworms * mhapc, :]
        elif fnum_recomb[i] == 0:
            fhapc = int(rand()/RAND_MAX)
            print(str(fhapc))
            h2[i, :] = fem[iix_fem[i] + fnworms * fhapc, :]
        else: 
            pass
    return(h1)


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
    #cdef np.ndarray fec 
    cdef float rr
    #cdef int[:] mate_array = np.empty(np.sum(females))
    for host in hosts:
        #chost = dfAdult.meta.index[dfAdult.meta.hostidx == host].values
        ahost = dfAdult.meta.hostidx == host
        females = np.logical_and(ahost, dfAdult.meta.sex == 'F').values
        males = np.logical_and(ahost, dfAdult.meta.sex == 'M').values
        if np.sum(males) == 0 or np.sum(females) == 0:
            print('Either there are 0 males in host or zero females in host')
            continue
        else:
            fec = dfAdult.meta.fec[females].values
            mate_array = np.random.randint(0, np.sum(males), np.sum(females),
                    dtype=np.int)
            # Parallelize this
            for loc in range(locus):
                rr = recombination_rate[loc]
                if rr == 0:
                    pass
                else:
                    cfemales = np.vstack((dfAdult.h1[str(loc)][females, :],
                        dfAdult.h2[str(loc)][females, :]))
                    cmales = np.vstack((dfAdult.h1[str(loc)][males, :],
                        dfAdult.h2[str(loc)][males, :]))
                    mate_worms(mate_array, 
                            fec, 
                            dfAdult.pos[str(loc)],
                            basepairs[loc],
                            rr,
                            cfemales, 
                            cmales)
    return(dfAdult)
