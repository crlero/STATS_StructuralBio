Stats README file
=====================

The Stats Project is student-level project for Structural Biology Informatics and
Python Programming courses attended at the Universitat Pompeu Fabra (Barcelona).
The aim of this project is the calculus of the statistical potential for protein
structures, a key feature to assess the accuracy of protein models.

This README file is intended primarily for people interested in working
with protein models and needing to assess its accuracy.
This statools package is a stand-alone application [stats] and also includes modules with
functions, which can be imported into Python library.


Installation
=================

To build and install Stats, download and unzip the source code, go to this
directory at the command line, and type::

    python setup.py build
    sudo python setup.py install

You might substitute ‘python’ with your specific version, for example ‘python3’.

Python Requirements
===================

We currently recommend using Python 3.4 from http://www.python.org, which
is the final version of Python 3. 


statools is currently supported and tested on the following Python versions:

- Python 3.4


Dependencies
============

statools needs in general to have installed the following dependencies:

- Biopython, see http://biopython.org /
  This package is mainly used to parse PDB files and to handle all the 
  biological data in order to gratefully carry out the neighbour-joining 
  algorithm here provided.
  
- NumPy, see http://www.numpy.org/
  This package is only used in the computationally-oriented modules.
  It is required for Bio.PDB module of Biopython.

- Matplotlib, see http://matplotlib.org/
  This package is used to visualize and plot the results.


Installation
============

First, make sure that Python 3.4 is installed correctly, as well as  NumPy and 
Matplotlib (see above) Packages, before installing statools.

Installation from source should be as simple as going to the statools
source code directory, and typing::

    python setup.py build
    sudo python setup.py install

If necessary, substitute 'python' with your specific version, for example 'python3'.


How to use this tool
=========================

statools is run by two Bioinformatic students from Universitat Pompeu Fabra (Barcelona).
We are willing to improve our coding skills, and whatever else comes up in this hybrid
field of Informatics and Biolgy.

Python Package for the Calculus of Statistical Potentials for Protein
Structure. stats calculates the statistical potential for each residue of the
protein and outputs either the energy/residue or a nice-visual plot of
energies. Calculus are done making use of another Python third-packages and modules
(Biopython, numpy, MatPlotlab, etc.) and following a neighbour-joining
algorithm. 
Authors: Jorge Roel and Cristina Leal. Year: 2015.

optional arguments:
  -h, --help            show this help message and exit
  -i [INFILE], --input [INFILE]
                        Input file. Input file needs to have a PDB Format
                        containg protein atoms and its coordinates (either
                        *.pdb or *.ent format).
  -o [OUTFILE], --output [OUTFILE]
                        Outfiles the Position-Residue and its associated
                        statistical potential. If not defined, results are
                        printed to standard output.
  -r [RADIUS], --radius [RADIUS]
                        Defines the desired radius of neighbor-joining
                        algorithm. If not set, the default radius given is
                        10Å. Radius has to be set to either: 5, 10, 15, 20, or
                        25Å.
  -v, --visual          If this argument is defined, results are plotted and
                        saved as 'filename.png'. If not used, plot is not
                        showed

From the command line:
usage: stats [-h] [-i [INFILE]] [-o [OUTFILE]] [-r [RADIUS]] [-v]

If you have any query, please contact us at: jorge.roel01@estudiant.upf.edu or
cristina.leal01@estudiant.upf.edu

Distribution Structure
======================

- README    		 -- This file.
- setup.py   		 -- Installation file.
- get_normalized.py	 -- Script used to build up the normalized dictionaries among 					1.110 protein domains from the PDB with <40% homology (to 				avoid family structure redundancy).
- statools		 -- Package folder
- stats	   		 --  Executable script 
- functions.py		 -- Main functions used for the executable. Can also be imported 				to the python package-library
- plot.py		 -- Main function used for plotting and save visual results