import unittest

import pandas as pd

from figs.worm import Worm
from figs.mutation import mutation_fx

class Test_Mutation(unittest.TestCase):
    def setUp(self):
        meta = pd.DataFrame()

        worms = Worm(meta)

        worms, newpositions = mutation_fx(2, worms, [0, 0.01],
                [200, 200])

if __name__ == '__main__':
    unittest.main()
