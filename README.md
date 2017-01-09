<!-- README.md is generated from README.Rmd. Please edit that file -->
rospca
======

The rospca package contains the implementation of robust sparse PCA using the ROSPCA algorithm of Hubert et al. (2016). Moreover, the simulation study and glass dataset discussed in this paper are included.

This package relies heavily on the code from Valentin Todorov for *rrcov* and on the *mrfDepth* package.

The latest development version of *rospca* can be installed from GitHub using

    install.packages("devtools")

    devtools::install_github("TReynkens/rospca")

If you work on Windows, make sure first that [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is installed when installing the development version from GitHub.

References
----------

Hubert, Mia, Tom Reynkens, Eric Schmitt, and Tim Verdonck. 2016. “Sparse PCA for High-Dimensional Data With Outliers.” *Technometrics* 58 (4): 424–34.
