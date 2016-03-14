<!-- README.md is generated from README.Rmd. Please edit that file -->
rospca
======

The rospca package contains the implementation of Robust Sparse PCA using the ROSPCA algorithm of Hubert et al. (201x). Moreover, the simulation study and glass dataset discussed in this paper are included.

This package relies heavily on the code from Valentin Todorov for *rrcov* and on the *mrfDepth* package.

The package is not available on CRAN but you can install the package using the binaries for *rospca* and *mrfDepthLight* from <https://wis.kuleuven.be/stat/robust/software#rospca>. The latest development version can be installed from GitHub using

    install.packages("devtools")

    devtools::install_github("TReynkens/rospca")
