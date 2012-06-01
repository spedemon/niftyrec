
from distutils.core import setup
import sys

VERSION = '1.6'
print "Python wrapper version: ",VERSION


setup(
    name = "NiftyRec",
    packages = ["NiftyRec"],
    version = VERSION,
    description = "Tomographic Reconstruction Toolbox, PET, SPECT, CT",
    author = "Stefano Pedemonte",
    author_email = "s.pedemonte@cs.ucl.ac.ul",
    url = "http://niftyrec.sourceforge.net/",
    download_url = "http://sourceforge.net/projects/niftyrec/files/nifty_rec-1.6.tar.gz",
    keywords = ["PET", "SPECT", "emission tomography", "transmission tomography", "tomographic reconstruction"],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    long_description = """\
NiftyRec Tomography Toolbox
---------------------------
NiftyRec is developed at the Centre for Medical Image Computing, University College London.
It provides tools for Emission and Transmission Tomographic reconstruction. Projection, back-projection and core iterative reconstruction routines are written in the C programming language and computationally intensive functions have a GPU accelerated version based on NVidia CUDA. 
It includes MLEM, OSEM, MAPEM reconstruction algorithms and examples for SPECT and PET, tough it can 
be used for diverse tomograhpic reconstruction problems. """
)


