\name{Glass}
\alias{Glass}

\title{
Glass data
}

\description{
Glass data of Lemberge et al. (2000) containing  Electron Probe X-ray Microanalysis (EPXMA) intensities for different wavelengths of 16--17th century archaeological glass vessels.
This dataset was also used in Hubert et al. (2005).
}

\format{
A data frame with 180 observations and 750 variables. These variables correspond to EPXMA intensities for different wavelengths and are indicated by \code{V1}, \code{V2}, ..., \code{V750}.
}


\usage{
data(Glass)
}



\source{
Lemberge, P., De Raedt, I., Janssens, K. H., Wei, F., and Van Espen, P. J. (2000), ``Quantitative Z-Analysis of the 16--17th Century Archaelogical Glass Vessels using PLS Regression of EPXMA and \eqn{\mu}-XRF Data," \emph{Journal of Chemometrics}, 14, 751--763.
}

\references{
Hubert, M., Rousseeuw, P. J., and Vanden Branden, K. (2005), ``ROBPCA: A New Approach to Robust Principal Component Analysis,'' \emph{Technometrics}, 47, 64--79.
}


\examples{
data(Glass)

res <- robpca(Glass, k=4, alpha=0.5)
matplot(res$loadings, type="l", lty=1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets  }% __ONLY ONE__ keyword per line