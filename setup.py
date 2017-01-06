import glob
from setuptools import setup

import numpy
import pysam


name = 'scottsmallsim'
version = '0.1'


try: 
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
    print("Cython not found")
else: use_cython = True


'''
if use_cython:
    print('******using cython******')
    ext_modules += [
            Extension("genda.transcripts.exon_utils",
                [""],
            ]
'''

metadata = {'name': name,
            'version': version,
            #'cmdclass': cmdclass,
            'scripts': glob.glob('scripts/*.py'),
            'author': 'Scott Smalls',
            }


if __name__ == '__main__':
    dist = setup(**metadata)
