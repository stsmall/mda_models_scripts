import unittest

import sys
import os.path
import random
import numpy as np
import pandas as pd
from pstats import Stats
import cProfile

from figs.recombination import recombination_locus, recombination_fx
from figs.worm import Worms
np.random.seed(20)


class Test_Recombination_Locus(unittest.TestCase):
    """Testing recombination_locus output
    """
    def setUp(self):
        self.basepairs = 10
        self.h1 = [1, 5, 6]
        self.h2 =  [2, 8]
        self.crossover_pos = 4

    '''
    def test_recombination_norm(self):
        h1, h2 = recombination_locus(self.h1,
                self.h2, self.crossover_pos)
        np.testing.assert_equal(h1, [1, 8])
        np.testing.assert_equal(h2, [2, 5, 6])

    def test_end(self):
        """Test when recombination happens at the end
        """
        h1, h2 = recombination_locus(self.h1,
                self.h2, 9)
        np.testing.assert_equal(h1, self.h1)
        np.testing.assert_equal(h2, self.h2)

    def test_start(self):
        h1, h2 = recombination_locus(self.h1,
                self.h2, 0)
        np.testing.assert_equal(h1, self.h2)
        np.testing.assert_equal(h2, self.h1)


    def test_multiple_crossovers(self):
        h1, h2 = recombination_locus(self.h1,
                self.h2, 4)
        h1, h2 = recombination_locus(h1, h2, 7)
        np.testing.assert_equal(self.h1, [1,5,6])
    '''


class Test_Recombination_Fx(unittest.TestCase):
    def setUp(self):
        villages = [0, 0, 0]
        sex = ['M', 'F', 'F']
        hostidx = ['v0h1', 'v0h1', 'v0h1']
        R0net = [0.5299222, 0.658231, 0.444]
        fec = [0, 10, 2]
        positions = {
                '0' : np.array([20, 30], dtype=np.uint64),
                '1' : np.array([1, 10, 50, 100], dtype=np.uint64)
                }

        loc0 = np.array([[0, 0],
                         [1, 1], 
                         [1, 0]],
                         dtype=np.uint8)

        hap1 = np.array([[0, 1, 0, 1], 
                         [0, 0, 1, 0], 
                         [1, 0, 0 , 0]], 
                         dtype = np.uint8)
        hap2 = np.array([[0, 0, 0, 0], 
                         [0, 1, 0, 1], 
                         [1, 1, 0 , 0]], 
                         dtype = np.uint8)

        meta = pd.DataFrame({
            'village': villages, 
            'sex': sex, 
            'hostidx': hostidx, 
            'fec': fec})
        worms = Worms(meta, haplotype1={'0' : loc0, '1' : hap1},
            haplotype2={'1': hap2}, positions=positions)
        self.worms = worms
        self.pr = cProfile.Profile()
        self.pr.enable()


    def tearDown(self):
        p = Stats(self.pr)
        p.strip_dirs()
        p.sort_stats('cumtime')
        #p.print_stats()

        

    def test_recombination_fx(self):
        df_adult_mf = recombination_fx(2, self.worms, [0, 0.5], [100, 200])
        #np.testing.assert_equal(self.adult.locus_0_h1[0] , [1, 3, 9])
        #np.testing.assert_equal(self.adult.locus_0_h2[0] , [2, 4, 6])
        #np.testing.assert_equal(df_adult_mf.shape[0], 12)






if __name__ == '__main__':
    unittest.main()
