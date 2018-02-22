\name{correctBatch}
\alias{correctBatch}
\title{
Batch correction approach based on the PPCCA model.
}
\description{
Batch effect is corrected by subtracting its effect only in the affected dimensions or PCs. 
Both the corrected and uncorrected dimensions are then used to predict the data.
}
\usage{
correctBatch(res11, designX, Ys, batchCL, grps, conF, nt, cGsub, comres)
}
\arguments{
  \item{res11}{
  Results of fitting the PPCCA model. 
}
  \item{designX}{
  Design matrix for the covariates
}
  \item{Ys}{
 Expression data
}
  \item{batchCL}{
  A vector identifying batches.
}
  \item{grps}{
  Number of batches
}
  \item{conF}{
  A vector identifying biological classes.
}
  \item{nt}{
  Number of biological classes
}
  \item{cGsub}{
  Colors
}
  \item{comres}{
  Results of ComBat correction
}
}
\details{
Batch effect is corrected by subtracting its effect only in the affected dimensions or PCs. 
Both the corrected and uncorrected dimensions are then used to predict the data.
}
\value{
  Results in correctBatch folder in the working directory.
}
\references{
Nyamundanda, G. and Sadanandam, A., 2016. A novel and robust statistical method to diagnose and correct batch effects in large-scale genomic data.
}
\author{
Gift Nyamundanda and Anguraj Sadanandam
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line