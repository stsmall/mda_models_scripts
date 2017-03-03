import unittest

import numpy as np
import pandas as pd
from pstats import Stats
import cProfile

from figs.recombination import recombination_fx
from figs.worm import Worms
from IPython import embed
np.random.seed(20)


class Test_Recombination_Locus(unittest.TestCase):
    """Testing recombination_locus output
    """
    def setUp(self):
        self.mate_array = np.array([
            0, 1], dtype=np.int64)
        self.fec = np.array([3], dtype=np.int64)



class Test_Recombination_Fx(unittest.TestCase):
    def setUp(self):
        worms = pd.read_pickle('')


    def tearDown(self):
        p = Stats(self.pr)
        p.strip_dirs()
        p.sort_stats('cumtime')
        #p.print_stats()
        

    def test_recombination_fx(self):
        # Female (index 3) mates with male 2 (index 1)
        df_adult_mf = recombination_fx(2, self.worms, [0, 0.005], [100, 200])
        self.assertEqual(df_adult_mf.meta.shape[0] , 12)
        self.assertEqual(df_adult_mf.h1['0'].shape[0], 12)
        self.assertEqual(df_adult_mf.h2['1'].shape[0], 12)
        np.testing.assert_equal(df_adult_mf.meta.index.values,
                                np.arange(12))
        # Recombination happens at index 3, 4, 4 (but 4 is at end so no
        # changes)
        np.testing.assert_equal(df_adult_mf.h1['1'][0, :] , [4, 4, 4, 0])

        # Test that the input array is not altered
        self.assertEqual(self.worms.h1['1'].shape[1], 4)

    def test_recombination_mutation(self):
        # Make sure that mutation doesn't alter the new worms.object
        pass

    





if __name__ == '__main__':
    unittest.main()
