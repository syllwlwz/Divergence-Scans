import numpy
import os
from distutils.core import setup,Extension

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

ext_modules=[ Extension('fastphase.fastphaseCython',["fastphase/fastphaseCython.c"],include_dirs=[numpy.get_include()]),
              Extension('fastphase.fastphaseCythonMT',["fastphase/fastphaseCythonMT.c"],include_dirs=[numpy.get_include()])]
setup(
    name='hapflk',
    version='1.3.0',
    description='haplotype-based test for differentiation in multiple populations',
    long_description=read('README.txt'),
    license = "GPL v3",
    author='Bertrand Servin',
    author_email="bertrand.servin@toulouse.inra.fr",
    url='https://forge-dga.jouy.inra.fr/projects/hapflk',
    packages=['fastphase','hapflk'],
    ext_modules = ext_modules,
    scripts=["hapflk/hapflk"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        ]
)
