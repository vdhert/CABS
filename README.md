![cabsdock-logo.jpg](https://bitbucket.org/repo/Bgk4jed/images/288308575-cabsdock-logo.jpg)
# README #
This is a repository for CABS-dock standalone application for protein-peptide molecular docking. CABS-dock enables protein-peptide docking with significant flexibility of protein receptor structure.


### Literature and general information on CABS-dock method ###

CABS-dock method was made available as a web server at [http://biocomp.chem.uw.edu.pl/CABSdock/](http://biocomp.chem.uw.edu.pl/CABSdock/). The following papers describe the CABS-dock web server and its example applications:

* CABS-dock web server for flexible docking of peptides to proteins without prior knowledge of the binding site, Nucleic Acids Research, 43(W1): W419-W424, 2015
* Modeling of protein-peptide interactions using the CABS-dock web server for binding site search and flexible docking, Methods, 93, 72-83, 2016
* Protein-peptide molecular docking with large-scale conformational changes: the p53-MDM2 interaction, Scientific Reports 6, 37532, 2016
* Highly flexible protein-peptide docking using CABS-dock, Methods in Molecular Biology, 1561: 69-94, 2017

* [Methods 93:72-83, 2016](http://www.sciencedirect.com/science/article/pii/S1046202315300207) - includes CABS-dock web server tutorial on visualization of docking results using VMD
* [Scientific Reports, 6:37532, 2016](https://www.nature.com/articles/srep37532) - describes modeling with large structural changes of protein receptor structure.
* [Methods in Molecular Biology, 1561:69-94, 2017](https://link.springer.com/protocol/10.1007%2F978-1-4939-6798-8_6) - includes tutorial on contact maps using CABS-dock web server

CABS-dock pipeline consist of the three following modules:

* Simulation module – performs docking simulations using peptide sequence, protein structure and set of parameters as an input. With default settings the module outputs a set of 10’000 of models (10 trajectories consisting of 1000 models) in C-alpha representation.

* Scoring module – selects representative and best-scored models from the simulation module output. Scoring module outputs sets of 10, 100 and 1000 top-scored models in C-alpha representation.

* Reconstruction to all-atom representation module – uses a Modeller package to reconstruct a set of 10 top-scored models from C-alpha to all-atom representation.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact