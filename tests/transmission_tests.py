import unittest

from IPython import embed

import sys
import os.path
sys.path.append(
            os.path.abspath(os.path.join(os.path.dirname(__file__),
                os.path.pardir)))

from transmission import transmission_fx


def main():
    """
    """
    dfMF = pd.DataFrame()
    dfJuv = pd.DataFrame() 
    hm = transmission_fx(2, [200, 200], 50, [20, 20], [8, 8], [True, True], [False, False], 
                      [12, 12], [36, 36], 1, dfMF, dfJuv, dfHost, deathdict) 
    embed()


if __name__ == '__main__':
    main()
