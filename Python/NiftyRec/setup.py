
from distutils.core import setup
import sys

try:
  VERSION_MAJ = '1'
  VERSION_MIN = '5'
  VERSION = str(VERSION_MAJ)+'.'+str(VERSION_MIN)
except:
  VERSION = 'x.x'

print "Python wrapper version: ",VERSION

setup(name='NiftyRec',
      version= VERSION,
      description='NiftyRec Python API',
      author='Stefano Pedemonte',
      author_email='s.pedemonte@cs.ucl.ac.uk',
      url='http://niftyrec.sourceforge.net/',
      py_modules = ['NiftyRec'],
     )



