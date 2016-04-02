Whole-cell modeling tutorial

== Contents
- About
- Requirements
- Installation
- Usage
- Further information
- Questions & feedback

== About
This tutorial is designed to teach students several of the fundamental concepts of whole-cell (WC) modeling including
- Model building,
- Multi-algorithm simulation,
- Parameter estimation, and
- Model annotation.

== Required background
The following background knowledge is needed to complete this tutorial:
- Flux balance analysis
- Stochastic simulation algorithm
- Python programming

== Required software
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

== Installation
To install the tutorial, simply download this repository to your computer. The following command can be used to retrieve the repository from GitHub:
```
git clone .git

```

== Using the tutorial
We recommmend
1. Watch the [introductory lecture](). This introduces WC modeling and provides an overview of the current state of the art. The [slides](Introduction.pdf) for this video are included in this repository.
2. Watch the [tutorial introduction video](). This breifly reviews WC modeling and outlines the tutorials.
3. Complete the tutorials 1-5. We recommend completing the tutorials in this order because the tutorials build upon each other. The [tutorial introduction slides] contain instructions for the tutorial and summarize the solutions. All of the materials for the tutorials, as well as their detailed solutions are contained in the directories titled "Exercise X ...".

== Further information
Please see these research articles and reviews for more information about WC modeling:
- Carrera Trends
- Karr Curr Opin Microbiol
- Karr et al. Cell
- Macklin Curr Opin Biotechnol
- Tomita

== Questions & feedback
This tutorial was developed by [Jonathan Karr](http://www.karrlab.org). Please contact [karr@mssm.edu](mailto:karr@mssm.edu) with any questions or feedback.
