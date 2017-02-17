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
from cpython cimport array
from libc.stdlib cimport rand, RAND_MAX
#from cython.parallel import parallel, prange
#from libc.stdlib cimport abort, malloc, free

DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t

#@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef np.ndarray[dtype=np.uint8_t, ndim=2] mate_worms(
        long[:] mate_array,
        np.ndarray[np.int64_t, ndim=1] fec,
        np.ndarray[np.int64_t, ndim=1] pos,
        int basepairs,
        float recomb_rate,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] fem,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] males):
    """ Mates and recombines at a given loci
    """
    # :TODO need to check max integer
    cdef np.intp_t i, j, l 
    cdef np.int64_t outsize
    cdef int k 
    outsize = np.sum(fec) 
    cdef int mnworms, fnworms
    cdef int mhapc, fhapc, recomb_pos 
    cdef long[:] iix_ma = np.repeat(mate_array, fec)
    cdef long[:] femindex = np.arange(fem.shape[0]/2, dtype=np.int64)
    cdef np.ndarray iix_fem = np.repeat(femindex, fec)
    cdef np.ndarray mnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    cdef np.ndarray fnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    cdef DTYPE_t[:, ::1] hout = np.empty((2*outsize, fem.shape[1]), dtype=np.uint8)
    mnworms = males.shape[0]/2
    fnworms = fem.shape[0]/2
    for i in range(outsize):
        print(mnum_recomb[i])
        if mnum_recomb[i] == 0:
            mhapc = int(rand()/RAND_MAX)
            hout[i, :] = males[iix_ma[i] + mnworms * mhapc, :]
        else:
            k = 0
            while k <= mnum_recomb[i]:
                k += 1
                recomb_pos = int(rand()/RAND_MAX*basepairs)
                pos_ix = 0
                for l in range(len(pos)):
                    if recomb_pos > pos[l]:
                        pos_ix = l
        if fnum_recomb[i] == 0:
            fhapc = int(rand()/RAND_MAX)
            hout[i + outsize, :] = fem[iix_fem[i] + fnworms * fhapc, :]
        else: 
            k = 0
            while k <= mnum_recomb[i]:
                k += 1
                recomb_pos = int(rand()/RAND_MAX*basepairs)
                for l in range(len(pos)):
                    if recomb_pos > pos[l]:
                        break
    return(hout)


def recombination_fx(locus,
                     dfAdult,
                     list recombination_rate,
                     list basepairs):
    """Calculate number of recombination events and rearranges haplotypes
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
    cdef np.ndarray out_array
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
                    dtype=np.int64)
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
                    out_array = mate_worms(mate_array,
                            fec,
                            dfAdult.pos[str(loc)],
                            basepairs[loc],
                            rr,
                            cfemales,
                            cmales)
    #new_meta = pd.DataFrame({})
    #dfAdult_mf = Worms(meta = pd
    return(out_array)
