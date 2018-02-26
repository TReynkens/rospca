\name{NEWS}
\title{News for package rospca}


\section{Changes in version 1.0.4 (2018-02-26)}{

 \itemize{
    \item Increase minimum required R version to 2.14.0 since byte-compilation of package code is only included since that version.
    \item Add ORCID ID.
 }
}

\section{Changes in version 1.0.3 (2017-10-20)}{

 \itemize{
    \item Enable byte-compilation of package code.
    \item Increase minimum required version of \pkg{mrfDepth} to 1.0.5.
 }
}

\section{Changes in version 1.0.2 (2017-04-03)}{

 \itemize{
    \item \code{adjOutlyingness} is renamed to \code{adjOutl} in \pkg{mrfDepth} version 1.0.3.
    The calls to this function are changed accordingly, and hence \pkg{mrfDepth} version 1.0.3 or higher is now required.
    \item Re-format R code.
    \item Fix bug in \code{robpca} when \code{k=1} and \code{skew=TRUE}.
 }
}


\section{Changes in version 1.0.1 (2017-02-28)}{

 \itemize{
    \item Change computation of score distances in skew-adjusted version of ROSPCA.
 }
}

\section{Changes in version 1.0.0 (2017-01-15)}{

 \itemize{
    \item Release on CRAN. The old version of \code{NEWS.Rd} can be found on GitHub: \url{https://github.com/TReynkens/rospca/blob/master/inst/NEWS.Rd}.
 }
}

