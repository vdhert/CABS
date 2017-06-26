# README #
This is a repository for CABS-dock standalone application for protein-peptide molecular docking. CABS-dock enables protein-peptide docking with significant flexibility of protein receptor structure.

**CABS-dock** pipeline consist of the three following modules:

* **Simulation module** – performs docking simulations using peptide sequence, protein structure and set of parameters as an input. With default settings the module outputs a set of 10’000 of models (10 trajectories consisting of 1000 models) in C-alpha representation.

* **Scoring module** – selects representative and best-scored models from a set 10’000 (in C-alpha representation) models. Outputs sets of 10, 100 and 1000 top-scored models

* **Reconstruction to all-atom representation module** – uses a Modeller package to reconstruct a set of top-scored models from C-alpha to all-atom representation

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