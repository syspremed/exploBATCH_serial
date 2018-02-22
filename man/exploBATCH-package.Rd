\name{exploBATCH-package}
\alias{exploBATCH-package}
\alias{exploBATCH}
\docType{package}
\title{
\packageTitle{exploBATCH}
}
\description{
\packageDescription{exploBATCH}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{exploBATCH}
\packageIndices{exploBATCH}
Explore batch effect (exploBATCH) is an R package for perfoming formal statistical tests to discover and correct batch effect. It searches for directions in the principal subspace that are associated with the batch variable of interest and correct for batch effect by subtracting this 
 effect only in the affected dimensions. Both the corrected and uncorrected dimensions are then used to predict the data.
}
\author{
\packageAuthor{exploBATCH}

Maintainer: \packageMaintainer{exploBATCH}
}
\references{
Nyamundanda, G. and Sadanandam, A., 2016. A novel and robust statistical method to diagnose and correct batch effects in large-scale genomic data.
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
## help(exploBATCH)
## Expression data and batch effect variable is required for this function.
data(Breast)
data(batchBreast)
expBATCH(D= Breast,batchBreast,Conf=NA,mindim=1,maxdim=10,method="ppcca")
}