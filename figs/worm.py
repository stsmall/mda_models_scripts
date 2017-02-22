"""
"""

from numpy import delete as ndelete
from numpy import vstack
import numpy as np
import pandas as pd


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


    def add_worms(self, df, index):
        """
        Parameters
        ----------
        df : figs.worm.Worms object
            other Worms object to add worms from
        index : int list
            numerical index from the other Worms object to add
        """
        try:
            self.meta = pd.concat([self.meta, df.meta.ix[index, :]], ignore_index=True)
            self.meta.reset_index(drop=True) #inplace=True
            # :TODO Make sure both self and df
            for i in df.h1.keys():
                try:
                    self.h1[i] = vstack((self.h1[i], df.h1[i][index,:]))
                except KeyError:
                    self.h1[i] = df.h1[i][index, :]
            for i in df.h2.keys():
                try:
                    self.h2[i] = vstack((self.h2[i], df.h2[i][index,:]))
                except KeyError:
                    self.h2[i] = df.h2[i][index, :]
        except ValueError:
            print("nothing to remove")

    def drop_worms(self, index):
        try:
            self.meta.drop(index)
            self.meta = self.meta.reset_index(drop=True) #inplace=True
            for i in self.h1.keys():
                self.h1[i] = ndelete(self.h1[i], index, axis=0)
            for i in self.h2.keys():
                self.h2[i] = ndelete(self.h2[i], index, axis=0)
        except ValueError:
            print("df empty")

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



