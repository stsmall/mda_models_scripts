import unittest

import sys
import os.path
sys.path.append(
            os.path.abspath(os.path.join(os.path.dirname(__file__),
                os.path.pardir)))

from recombination import recombination_locus, recombination_fx
import random
import numpy as np
import pandas as pd

random.seed(100)

class testRecombinationLocus(unittest.TestCase):
    """Testing recombination_locus output
    """
    def setUp(self):
        self.basepairs = 10
        self.h1 = [1, 5, 6]
        self.h2 =  [2, 8]
        self.crossover_pos = 4

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
        np.testing.assert_equal(h1, self.h1)
        np.testing.assert_equal(h2, self.h2)




class testRecombinationFx(unittest.TestCase):
    def setUp(self):
        villages = [0, 0]
        sex = ['M', 'F']
        hostidx = ['v0h1', 'v0h1']
        R0net = [0.5299222, 0.658231]
        fec = [0, 20]

        M_loc1 = [1, 3, 9]
        M_loc2 = []




if __name__ == '__main__':
    unittest.main()
