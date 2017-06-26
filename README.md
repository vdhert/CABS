# README #
This is a repository for CABS-dock standalone application for protein-peptide molecular docking. CABS-dock enables protein-peptide docking with significant flexibility of protein receptor structure.


### General information on CABS-dock method ###

CABS-dock method was made available as a web server in 2015 (see [Nucleic Acids Research 43:W419-W424, 2015](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkv456)) available at [http://biocomp.chem.uw.edu.pl/CABSdock/](http://biocomp.chem.uw.edu.pl/CABSdock/). Since that, the CABS-dock has been described in: 

* [Methods 93:72-83, 2016](http://www.sciencedirect.com/science/article/pii/S1046202315300207) - includes CABS-dock tutorial on visualization of docking results using VMD
* [Scientific Reports, 6:37532, 2016](https://www.nature.com/articles/srep37532) - describes modeling with large structural changes of protein receptor structure.
* [Methods in Molecular Biology, 1561:69-94, 2017](https://link.springer.com/protocol/10.1007%2F978-1-4939-6798-8_6) - includes tutorial on contact maps

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