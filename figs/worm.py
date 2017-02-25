"""
"""

from numpy import delete as ndelete
from numpy import vstack
import numpy as np
import pandas as pd



def merge_positions(pos1, pos2, newpos=None):
    """ Return indexes where to insert
    """
    # This could be sped up
    for i, j in enumerate(pos1):
        pass
    for i, j in enumerate(pos2):
        pass






class Worms(object):
    def __init__(self, meta, haplotype1=None, haplotype2=None,
            positions = None, selection=None, cds_coords=None):
        self.meta = meta
        if haplotype1:
            self.h1 = haplotype1
        else:
            self.h1= {}
        if haplotype2:
            self.h2 = haplotype2
        else:
            self.h2 = {}
        if positions:
            self.pos = positions
        else:
            self.pos = {}
        if selection:
            self.sel = selection
        else:
            self.sel = {}
        if cds_coords:
            self.coord = cds_coords
        else:
            self.coord = {}


    def _merge_positions(self, loc, oworm, newpos = None):
        # Not the fastest
        pos1 = self.pos[loc]
        pos2 = oworms.pos[loc]
        common = np.intersect1d(self.pos[loc], oworms.pos[loc])
        m1 = [i for i in pos1 if i not in common]
        m2 = [i for i in pos2 if i not in common]

        
        n1 = self.h1[loc].shape[0]
        for i in m1:
            iix = np.argmax(pos1 > i)
            np.insert(self.h1[loc], iix, 
                    np.zeroes(n1, dtype=np.uint8))
            np.insert(self.h2[loc], iix, 
                    np.zeroes(n1, dtype=np.uint8))
            positions = np.insert(pos1, iix, i)

        n2 = self.h1[loc].shape[0]
        for i in m2:
            iix = np.argmax(pos2 > i)
            np.insert(oworm.h1[loc], iix, 
                    np.zeroes(n2, dtype=np.uint8))
            np.insert(oworm.h2[loc], iix, 
                    np.zeroes(n2, dtype=np.uint8))
            positions = np.insert(pos2, iix, i)

        return(oworm)

            

    def add_worms(self, df, index):
        """
        Parameters
        ----------
        df : figs.worm.Worms object
            other Worms object to add worms from
        index : int list
            numerical index from the other Worms object to add
        """
        if len(index) != 0 and self.meat.shape[0] !=0:
            self.meta = pd.concat([self.meta, df.meta.ix[index, :]], ignore_index=True)
            self.meta.reset_index(drop=True) #inplace=True
            for i in df.h1.keys():
                try:
                    assert self.h1[i].shape[1] == df.h1[i].shape[1]
                    self.h1[i] = vstack((self.h1[i], df.h1[i][index,:]))
                except KeyError:
                    self.h1[i] = df.h1[i][index, :]
                    self.pos[i] = df.pos[i]
            for i in df.h2.keys():
                try:
                    assert self.h2[i].shape[1] == df.h2[i].shape[1]
                    self.h2[i] = vstack((self.h2[i], df.h2[i][index,:]))
                except KeyError:
                    self.h2[i] = df.h2[i][index, :]
        elif self.meta.shape[0] == 0 and len(index) != 0:
            for i in df.h1.keys():

        else:
            self.meta = pd.concat([self.meta, df.meta.ix[index, :]], ignore_index=True)
            self.meta.reset_index(drop=True) #inplace=True
            print("Nothing to add")


    def drop_worms(self, index):
        try:
            self.meta.drop(index)
            self.meta = self.meta.reset_index(drop=True) #inplace=True
            for i in self.h1.keys():
                self.h1[i] = ndelete(self.h1[i], index, axis=0)
            for i in self.h2.keys():
                self.h2[i] = ndelete(self.h2[i], index, axis=0)
        except ValueError:
            print("DF empty in drop")


    def calc_allele_frequencies(self, host=None, village=None):
        """ Calculate allele frequencies
        """
        all_loci_shape = [i.shape[0] for _, i in self.h1.iteritems()]
        allele_freqs = np.empty(np.sum(all_loci_shape), dtype=np.float64)
        for loc in self.h1.keys():
            if loc in self.h2.keys():
                allele_freq[loc] = np.sum(self.h1[loc] +\
                        self.h2[loc])/float(self.h1[loc].shape[0])
            else:
                allele_freq[loc] = np.sum(self.h1[loc])
