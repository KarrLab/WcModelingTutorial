# Whole-cell modeling tutorial

## Contents
- About
- Requirements
- Installation
- Usage
- Further information
- Questions & feedback

## About
This tutorial is designed to teach students several of the fundamental concepts of whole-cell (WC) modeling including
- Model building,
- Multi-algorithm simulation,
- Parameter estimation, and
- Model annotation.

## Required background
The following background knowledge is needed to complete this tutorial:

- Flux balance analysis
- Stochastic simulation algorithm
- Python programming

## Required software
The tutorial requires the following packages:

- biopython
- cobra
- libglpk-dev
- matplotlib
- numpy
- openpyxl
- python 2.7
- scipy

On Ubuntu, these can be installed with the following commands:
```
sudo apt-get install python numpy scipy matplotlib libglpk-dev biopython
sudo pip install openpyxl cobra
```

## Installation
To install the tutorial, simply download this repository to your computer. The following command can be used to retrieve the repository from GitHub:
```
git clone https://github.com/KarrLab/WcModelingTutorial.git

```

## Using the tutorial
We recommmend:

1. Watch the [introductory lecture](). This introduces WC modeling and provides an overview of the current state of the art. The [slides](https://github.com/KarrLab/WcModelingTutorial/raw/master/1.%20Introduction%20to%20whole-cell%20modeling.pdf) for this video are included in this repository.
2. Watch the [tutorial introduction video](). This breifly reviews WC modeling and outlines the tutorials. The [slides](https://github.com/KarrLab/WcModelingTutorial/raw/master/3.%20Exercises.pdf) for this video are included in this repository.
3. Complete the tutorials 1-5. We recommend completing the tutorials in this order because the tutorials build upon each other. The [tutorial introduction slides](https://github.com/KarrLab/WcModelingTutorial/raw/master/3.%20Exercises.pdf) contain instructions for the tutorial and summarize the solutions. All of the materials for the tutorials, as well as their detailed solutions are contained in the directories titled "Exercise X ...".

## Further information
Please see these research articles and reviews for more information about WC modeling:

- Carrera J & Covert MW. Why Build Whole-Cell Models?. *Trends Cell Biol* 25, 719–22 (2015). doi: [10.1016/j.tcb.2015.09.004](http://dx.doi.org/10.1016/j.tcb.2015.09.004)
- Karr JR, Takahasi K & Funahashi A. The principles of whole-cell modeling. *Curr Opin Microbiol* 27, 18&ndash;24 (2015). doi: [10.1016/j.mib.2015.06.004](http://dx.doi.org/10.1016/j.mib.2015.06.004)
- Karr JR, Sanghvi JC, Macklin DN, Gutschow MV, Jacobs JM, Bolival B Jr, Assad-Garcia N, Glass JI & Covert MW. A whole-cell computational model predicts phenotype from genotype. *Cell* 150, 389&ndash;401 (2012). doi: [10.1016/j.cell.2012.05.044](http://dx.doi.org/10.1016/j.cell.2012.05.044)
- Macklin DN, Ruggero NA & Covert MW. The future of whole-cell modeling. *Curr Opin Biotechnol* 28, 111&ndash;115 (2014). doi: [10.1016/j.copbio.2014.01.012](http://dx.doi.org/10.1016/j.copbio.2014.01.012)
- Tomita M, Hashimoto K, Takahashi K, Shimizu TS, Matsuzaki Y, Miyoshi F, Saito K, Tanida S, Yugi K, Venter JC et al. E-CELL: software environment for whole-cell simulation. *Bioinformatics* 15, 72&ndash;84 (1999). doi: [10.1093/bioinformatics/15.1.72](http://dx.doi.org/10.1093/bioinformatics/15.1.72)
- Tomita M. Whole-cell simulation: a grand challenge of the 21st century. *Trends Biotechnol* 19, 205&ndash;210 (2001). doi: [10.1016/S0167-7799(01)01636-5](http://dx.doi.org/10.1016/S0167-7799(01)01636-5)

## Questions & feedback
This tutorial was developed by [Jonathan Karr](http://www.karrlab.org). Please contact [karr@mssm.edu](mailto:karr@mssm.edu) with any questions or feedback.
