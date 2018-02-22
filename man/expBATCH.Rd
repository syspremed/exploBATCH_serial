\name{expBATCH}
\alias{expBATCH}
\title{
Discovering and correcting batch effect.
}
\description{
A tool for discovering and correcting for batch effect using an approach in Nyamundanda et al (2016), which is based on the PPCCA model.
}
\usage{
expBATCH(D, batchCL, Conf = NA, mindim = 1, maxdim = 2, method = "ppcca")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
Expression dataset in samples by features format.
}
  \item{batchCL}{
A vector identifying batches.
}
  \item{Conf}{
A vector identifying biological classes.
}
  \item{mindim}{
Minimum number of principal components to be considered. Default is 1.
}
  \item{maxdim}{
Maximum number of principal components to be considered. Default is 2. maxdim should be greater or equal to mindim.
}
  \item{method}{
Method for batch correction. Either PPCCA or ComBat. PPCCA is the default.
}
}
\details{
It searches for directions in the principal subspace that are associated with the batch effect and correct for batch effect by subtracting this 
effect only in the affected dimensions. Both the corrected and uncorrected dimensions are then used to predict the data.
}
\value{
The output of this pipeline is saved in your working directory.
}
\references{
Nyamundanda, G. and Sadanandam, A., 2016. A novel and robust statistical method to diagnose and correct batch effects in large-scale genomic data.
}
\author{
Gift Nyamundanda and Anguraj Sadanandam
}
\examples{
## help(expBATCH)
## Expression data and batch effect variable is required for this function.
data(Breast)
data(batchBreast)
expBATCH(D= Breast,batchBreast,Conf=NA,mindim=1,maxdim=10,method="ppcca")
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line