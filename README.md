# README
This is a repository for CABSdock standalone application for molecular docking of peptides to proteins. The CABSdock allows for flexible docking (also with large-scale conformational changes) without knowledge about the binding site. The CABSdock enables peptide docking using only information about the peptide sequence and the protein receptor structure. In addition to this input data, many advanced options are available that allow for manipulation of a simulation setup, a degree of protein flexibility or guiding the peptide binding etc.
###Detailed instructions and tutorials are provided on [CABSdock WIKI PAGE](https://bitbucket.org/lcbio/cabsdock/wiki)
### To install CABS on Linux Debian / Ubuntu / Mint:
####1. Install python-pip
```
#!bash
sudo apt install python-pip
```
####2. Download latest version of CABS
```
#!bash
wget https://bitbucket.org/lcbio/cabsdock/downloads/CABS-<version>.tar.gz
```
####3. Install CABS with pip (make sure you are using pip for python 2.7.*)
#####For all users
```
#!bash
sudo pip install CABS-<version>.tar.gz
```
#####Locally
```
#!bash
pip install CABS-<version>.tar.gz --user  
```
####4. To uninstall CABS with pip run:
#####For all users
```
#!bash
sudo -H pip uninstall CABS
```
#####Locally 
```
#!bash
pip uninstall CABS  
```
###To run CABSdock simply type:
```
#!bash
cabsDock
```
###Additionally, running CABSdock requires installation of the following packages:
#####Fortran compiler
(preferably GNU gfortran)

Install with:
```
#!bash
sudo apt install gfortran
```
For other compilers remember to run CABSdock with **"--fortran-compiler <your-compiler>"** option

#####DSSP
(optional if CABSdock will be run offline)

Install with:
```
#!bash
sudo apt install dssp
```
#####Modeller
(optional for all-atom reconstruction of models)

Follow [instructions](https://salilab.org/modeller/download_installation.html) on how to install it on your machine.

### Windows and Mac versions coming soon !!!

------------------------------------------------------------------------------------------------------------------------

# ABOUT THE METHOD #

CABSdock method has been first made available as a web [server](http://biocomp.chem.uw.edu.pl/CABSdock). The standalone application version [submitted to publication] provides the same modeling methodology equipped with many additional features and customizable options.

The following papers describe the CABS-dock method/ web server/ and its example applications:

* [CABS-dock web server for flexible docking of peptides to proteins without prior knowledge of the binding site, Nucleic Acids Research, 43(W1): W419-W424, 2015](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkv456)
* [Modeling of protein-peptide interactions using the CABS-dock web server for binding site search and flexible docking, Methods, 93, 72-83, 2016](http://www.sciencedirect.com/science/article/pii/S1046202315300207)
* [Protein-peptide molecular docking with large-scale conformational changes: the p53-MDM2 interaction, Scientific Reports 6, 37532, 2016](https://www.nature.com/articles/srep37532)
* [Highly flexible protein-peptide docking using CABS-dock, Methods in Molecular Biology, 1561: 69-94, 2017](https://link.springer.com/protocol/10.1007%2F978-1-4939-6798-8_6)

CABS-dock pipeline consist of the three following modules:

* Simulation module – performs docking simulations using peptide sequence, protein structure and set of parameters as an input. With default settings the module outputs a set of 10’000 of models (10 trajectories consisting of 1000 models) in C-alpha representation.
* Scoring module – selects representative and best-scored models from the simulation module output. Scoring module outputs sets of 10, 100 and 1000 top-scored models in C-alpha representation.
* Reconstruction to all-atom representation module – uses a Modeller package to reconstruct a set of 10 top-scored models from C-alpha to all-atom representation.

CABS-dock application uses some external software packages for the following tasks:

* [Gfortran](https://gcc.gnu.org/wiki/GFortran) - complier for the compilation of CABS simulation model code
* [dssp](http://www.cmbi.ru.nl/dssp.html) - program for secondary structure assignment of protein receptors from PDB files
* [MODELLER](https://salilab.org/modeller/) - program for modeling of protein structure using spatial restraints. CABSdock uses MODELLER for reconstruction of predicted complexes from C-alpha to all-atom representation. 
* [Matplotlib](https://matplotlib.org/) - Python 2D plotting library which produces publication quality figures. CABSdock uses Matplotlib for automated analysis of simulation results - several kinds of plots are being made. 
* [numpy](http://www.cmbi.ru.nl/dssp.html) - package for scientific computing with Python

------------------------------------------
Laboratory of Computational Biology, 2017