
from distutils.core import setup
import sys

try:
  VERSION_MAJ = '@Nifty_Rec_VERSION_MAJOR@'
  VERSION_MIN = '@Nifty_Rec_VERSION_MINOR@'
  VERSION = str(VERSION_MAJ)+'.'+str(VERSION_MIN)
except:
  VERSION = ''

print "Python wrapper version: ",VERSION

setup(name='NiftyRec',
      version= VERSION,
      description='NiftyRec Python API',
      author='Stefano Pedemonte',
      author_email='s.pedemonte@cs.ucl.ac.uk',
      url='http://niftyrec.sourceforge.net/',
      py_modules = ['NiftyRec.NiftyRec'],
     )



