# Complete Moveset Optimization Algorithm

This repository contains Python scripts for the Complete Moveset Search (CMS) optimization algorithm. This algorithm can be used for general optimization of single objective functions.
This algorithm was developed for fitting a biophysical and population genetic model for the analysis of the codon usage bias described in [this](https://www.biorxiv.org/content/10.1101/578815v1) paper.
A detailed description of the algorithm can be found in the Supplemental Material.

If you use this algorithm for a project or published research, please cite 
##### Kion-Crosby WB, Manhart M, Morozov AV. 2019. Inferring biophysical models of evolution from genome-wide patterns of codon usage. bioRxiv doi: 10.1101/578815 

# Installation

To install the optimizer, first clone the repository into a desired location. Move to the `Complete_Moveset_Search` directory in the terminal and run

`pip install -e .`

This will install the CMS optimization algorithm. Once the installation is complete, the algorithm can be accessed by Python in any directory.

# Running the Optimizer

Once installed, the algorithm can be run to find the optimum (minimum) of a standard Python function after import:

`from CMS import optimize`

An example of this can be found in the `Examples` folder.

