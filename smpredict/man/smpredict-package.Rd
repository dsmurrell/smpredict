\name{smpredict-package}
\alias{smpredict-package}
\alias{smpredict}
\docType{package}
\title{
Predicts various properties on small drug-like molecules 
~~ smpredict ~~
}
\description{
%% More about what it does (maybe more than one line)
%% ~~ A concise (1-5 lines) description of the package ~~

smpredict is a package that, given a set of
small molecules, predicts various important properties of these
molecules. At the moment the only property type available is logP -
the log of the ratio of concentrations of a compound in the two phases
of a mixture of two immiscible solvents at equilibrium.
}
\details{
\tabular{ll}{
Package: \tab smpredict\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2011-12-08\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Daniel Murrell <dsmurrell@gmail.com>
}
\references{
TBD: insert reference to published paper
TBD: insert reference to the caret publication
TBD: insert reference to the PD publication
}
\keyword{ logP }
\seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
predictions <- PredictLogP(system.file("test_structures", "structures_10.sdf", package="smpredict"))
}
