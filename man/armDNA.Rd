\name{armDNA}
\alias{armDNA}
\title{
Analysis of the effects of age and covariates on mDNA profiles.
}
\description{
This function estimates and tests the effects of age and covariates on mDNA profiles.
The transformation function for the age effect on the mDNA profiles can be either estimated nonparametrically or specified to be a power transformation function. Multple CPUs can be utilized through parallel computation.}
\usage{
armDNA(dat, pow=NULL, age.num=5, cl.cores=1)
}
\arguments{
  \item{dat}{
A list contains three variables, namely age (x), mDNA profile (y) and covariate (z), where x is a n-vector (the sample size), y is a  n x m matrix (m: the number of mDNA marks), z is a n x p matrix (p: the number of covariates).
}
  \item{pow}{
This can be used to specify the transformation function for the age effect, which should be either NULL (unknown piecewise linear transformation will be estimated) or a positive number (a power transformation will be specified).
}
  \item{age.num}{
The number of nodes for the age effects, which should be an integer greater than 2, with the default value 5. This parameter should be specified only if "pow" = NULL.
}
  \item{cl.cores}{
The number of CPU cores used in parallel computation, which should be a positive integer, with the default value 1 (see the help page for the R package "parallel").
}
}
\value{
\item{par}{Estimates of the effects of age and covariates on methylation profile. Column 1 is for the age effects, and columns 2 through p+1 are for covariate effects}
\item{se}{Standard errors of "par"}
\item{p.value}{P-values for Wald test of the effects of age and covariates on methylation profile}
\item{age.h}{Estimated transformation function result. If "pow" is a positive number, "age.h" is NULL; if "pow"=NULL, "age.h" is a age.num x 2 matrix, with the first column being the age values at age.num nodes and the second column being the corresponding estimated H values.}
}
\author{
Chenyang Wang and Hong Zhang
}

\examples{
    data(data);
    ###########
    # Not run:
    # res = armDNA(data);
}

\keyword{methylation}
\keyword{age}
\keyword{association}
\keyword{Wald test}
