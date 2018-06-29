 
These codes recreate the experiments of the NIPS paper "regularized nonlinear acceleration" (authors: Damien Scieur, Alexandre d'Aspremont, Francis Bach), using Matlab.

Link to the paper: 
- (direct download)     http://papers.nips.cc/paper/6267-regularized-nonlinear-acceleration.pdf
- (NIPS.cc)             http://papers.nips.cc/paper/6267-regularized-nonlinear-acceleration

Main contact: damien.scieur@inria.fr


There are two main folders:

- basic_example: Contains a standalone file "rmpe_logistic_basic.m", which shows how to plug RMPE algorithm to the gradient method, and compare the algorithm with AMPE (extrapolation without regularization).

- advanced: Contains a package for numerical experiments involving comparisons with RMPE. If you want to recreate the experiments of the paper, include "advanced" folder and all of its subfolders in the matlab path and run "solve_logistic.m".


BibTex:

@incollection{NIPS2016_6267,
title = {Regularized Nonlinear Acceleration},
author = {Scieur, Damien and d\textquotesingle Aspremont, Alexandre and Bach, Francis},
booktitle = {Advances in Neural Information Processing Systems 29},
editor = {D. D. Lee and M. Sugiyama and U. V. Luxburg and I. Guyon and R. Garnett},
pages = {712--720},
year = {2016},
publisher = {Curran Associates, Inc.},
url = {http://papers.nips.cc/paper/6267-regularized-nonlinear-acceleration.pdf}
}