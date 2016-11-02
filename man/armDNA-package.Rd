\name{armDNA-package}
\alias{armDNA-package}
\docType{package}
\title{
Age-related DNA methylation
}
\description{
This R package is developed to detect age-related genomewide DNA methylation marks based on a functional beta model using mDNA data obtained through the widely used bisulfite conversion technique. Let mDNA bs the beta value defined as (M+50)/(M+U+100), where M and U are methylated profile and unmethylated profile, respectively, then the mean mDNA profile is related to the age and covariates through the following functional beta model:

logit{E(mDNA)} = beta0 + beta1 * H(Age) + beta2 * Covar,

where H() is an increasing transformation function, and Covar is a covariate vector. The unknown parameters beta0, beta1, beta2 are mDNA-mark dependent. H() is independent of mDNA-marks, which can either be estimated nonparametrically or be specified to be a power function. 

H() is assumed to be a piecewise linear function with K nodes (K sample quantiles of ages) if H() is estimated; otherwise, H() is a prespecified power transformation function.

Multiple CPUs can be used to speed up the algorithm through parallel computation.

}
\details{
\tabular{ll}{
Package: \tab armDNA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-08-03\cr
License: \tab Artistic License 2.0\cr
}
}
\author{
Chenyang Wang and Hong Zhang

Maintainer: Hong Zhang <zhanghfd@fudan.edu.cn>
}
\references{
Wang C, Shen Q, Du L, Xu J, and Zhang H (2016). armDNA: A functional beta model for detecting age-related genomewide DNA methylation marks
}
\keyword{package}
