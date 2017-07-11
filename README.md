# README #
This is a repository for CABSdock standalone application for protein-peptide molecular docking. CABSdock method is also available as a web server at [http://biocomp.chem.uw.edu.pl/CABSdock/](http://biocomp.chem.uw.edu.pl/CABSdock/).
## Detailed instructions and tutorials are provided on [CABSdock WIKI PAGE](https://bitbucket.org/lcbio/cabsdock/wiki/) ##

-------------------------------------------

# INSTALLATION #
CABSdock standalone package ... full functionality requires installation of the following software packages:

* 
* [MODELLER](https://salilab.org/modeller/) - a program for comparative modeling of protein structure using spatial restraints. CABSdock uses MODELLER for reconstruction of predicted complexes from C-alpha to all-atom representation. 
* [Matplotlib](https://matplotlib.org/) - a Python 2D plotting library which produces publication quality figures. CABSdock uses Matplotlib for automated analysis of simulation results - several kinds of plots are being made. 
* ...

--------------------------------------------

# ABOUT THE METHOD ###

The following papers describe the CABS-dock method/ web server/ and its example applications:

* [CABS-dock web server for flexible docking of peptides to proteins without prior knowledge of the binding site, Nucleic Acids Research, 43(W1): W419-W424, 2015](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkv456)
* [Modeling of protein-peptide interactions using the CABS-dock web server for binding site search and flexible docking, Methods, 93, 72-83, 2016](http://www.sciencedirect.com/science/article/pii/S1046202315300207)
* [Protein-peptide molecular docking with large-scale conformational changes: the p53-MDM2 interaction, Scientific Reports 6, 37532, 2016](https://www.nature.com/articles/srep37532)
* [Highly flexible protein-peptide docking using CABS-dock, Methods in Molecular Biology, 1561: 69-94, 2017](https://link.springer.com/protocol/10.1007%2F978-1-4939-6798-8_6)

CABS-dock pipeline consist of the three following modules:

* Simulation module – performs docking simulations using peptide sequence, protein structure and set of parameters as an input. With default settings the module outputs a set of 10’000 of models (10 trajectories consisting of 1000 models) in C-alpha representation.

* Scoring module – selects representative and best-scored models from the simulation module output. Scoring module outputs sets of 10, 100 and 1000 top-scored models in C-alpha representation.

* Reconstruction to all-atom representation module – uses a Modeller package to reconstruct a set of 10 top-scored models from C-alpha to all-atom representation.


--------------------
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)