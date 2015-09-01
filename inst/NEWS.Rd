\name{NEWS}
\title{News for package rospca}

\section{Changes in rospca version 0.2-00 (2015-09-01)}{
  \itemize{
      \item Changes in R code:
      \itemize{
          \item Included new standardisation scheme and SD reweighting.
          \item Added option to compute cutoff for ODs in \code{distPCA}, default is \code{TRUE}.
          \item \code{H0} is now in the output list of \code{rospca}.
          \item Outlyingness part now checks if the obtained outlyingness is valid.
          \item Updated function calls to \code{outlyingness} and \code{adjOutlyingness} for compatibility with \pkg{mrfDepth} versions \eqn{\ge 0.3.05}. Still works for older versions of \pkg{mrfDepth}.
          \item Added error message when scales are 0 in \code{matStand}.
          \item Included warning in \code{rospca} when \code{ndir="all"} and there are more than 500 observations.
          \item Fixed bug with indexset in \code{robpca}.
          \item Set xlim and ylim in \code{diagPlot}.
          
      }
      
     \item Documentation changes:
      \itemize{
          \item Updated documentation of \code{rospca}.
          \item Updated documentation of \code{robpca}.
      }
      
    \item Misscelaneous changes:
      \itemize{
          \item Updated title to have correct upper case format.
          \item Import default packages as obliged by R 3.3.0.
          \item Do not import \pkg{compiler} anymore.
      }
  }
}

\section{Changes in rospca version 0.1-06 (2015-05-12)}{
  \itemize{
      \item Changes in R code:
      \itemize{
          \item \code{H0} is now in the output list of \code{rospca_part1} and is also
          an input argument of \code{rospca_part2} now.
          \item Changes with \code{H0} also incorporated in \code{selectLambda}.
      }
  }
}


\section{Changes in rospca version 0.1-05 (2015-05-06)}{
  \itemize{
      \item Changes in R code:
      \itemize{
          \item \code{selectLambda} now has an \code{ndir} option that is used in \code{rospca}.
          \item The default value for \code{kmax} in \code{selectLambda} is now 10. The argument \code{k} for \code{selectLamb
          da} does not have a default value anymore.
          \item Fixed problem with \code{classPC} in \code{robpca} (line 327).
      }
  }
}



\section{Changes in rospca version 0.1-04 (2015-03-24)}{
  \itemize{
      \item Changes in R code:
      \itemize{
          \item Fixed problem with \code{switch} in \code{select_method} (not exported).
      }
  \item Miscellaneous changes:
    \itemize{
    \item Changed NEWS file to R documentation style.
    }
  }
}


\section{Changes in rospca version 0.1-03 (2015-03-05)}{
  \itemize{
    \item Changes in R code:
      \itemize{
        \item Using \code{ndir="all"} in \code{ROSPCA} and \code{ROBPCA} when \code{n>500} now gives a warning.
    
        \item \code{H0} is now a logical vector in \code{ROBPCA} and \code{ROSPCA}.
    
	      \item Renamed \code{datagen} to \code{dataGen}.
     }
    \item Documentation changes:
      \itemize{
          \item Documentation update.
      }
   \item Miscellaneous changes:
    \itemize{
  	  \item Updated format of package title to R 3.2.0 rules.
   }
  }
}




\section{Changes in rospca version 0.1-02 (2015-02-11)}{
  \itemize{
  \item Documentation changes:
      \itemize{
      	\item Changed keywords.
      	
      	\item Updated 2014 to 2015 for Hubert et al. paper.
      	
      	\item Documentation update.
      }
  }
}


\section{Changes in rospca version 0.1-01 (2014-12-04)}{
  \itemize{
    \item Changes in R code:
      \itemize{
      	\item Updated \code{robpca} and \code{rospca} to use \code{classPC} from \pkg{robustbase}
  	        instead of \code{kernelSVD} and \code{kernelEVD} from \pkg{rrcov}.
  	  
      	\item Removed calls to \code{repmat} and \code{zeros} in \code{robpca}.
      	 
      	\item Included \code{measure} in output from \code{selectLambda}.
      	
      	\item Included univariate MCD function (\code{unimcd}, not exported).
      	
      	\item Changed code style when necessary.
     }
    \item Documentation changes:
      \itemize{
          \item Documentation update.
      }
   \item Miscellaneous changes:
    \itemize{
  	  \item Changed necessary \pkg{rrcovHD} version to 0.2-3.
   }
  
  }
}
	




  



 