\name{NEWS}
\title{News for package rospca}


\section{Changes in rospca version 0.5.0 (20xx-yy-zz)}{
  \itemize{
    \item Add variable names as row names for loadings matrix.
    \item Add row names for scores matrix.
    \item Add names for \code{od} and \code{sd}.
    \item Add names for eigenvalues and center.
    \item Add names for \code{H0}, \code{H1}, \code{H2} and \code{H3}.
    \item Add names for \code{D}.
    \item Remove byte-compiler.
    \item Reformat DESCRIPTION file.
  }
}

\section{Changes in rospca version 0.4.0 (2016-12-01)}{
  \itemize{
    \item Do not apply ROBPCA with MCD when \code{skew=TRUE}.
    \item Fix error with indexset in ROBPCA with MCD.
    \item Return fit corresponding to optimal \eqn{\lambda} in \code{selectLambda}.
    \item Set needed \pkg{robustbase} version to \emph{0.92-6} to solve problems with \code{k=1} in \code{robpca}.
    \item Update Hubert et al. (2016) reference.
    \item Improve documentation.
    \item Change version numbering format.
    \item Use \pkg{mrfDepth} instead of \pkg{mrfDepthLight}.
  }
}

\section{Changes in rospca version 0.3.0 (2016-03-14)}{
  \itemize{
    \item Release.
  }
}

