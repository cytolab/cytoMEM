\name{MEM_matrix}
\alias{MEM_matrix}
\docType{data}
\title{
MEM matrix
}
\description{
This matrix is the output generated from \code{\link{MEM}} analysis of the \code{\link{PBMC}} dataset. It is meant to be used as input for the \code{\link{MEM_RMSD}} function to generate RMSD scores of similarity.
}
\usage{data(MEM_matrix)}
\format{
  The format is the 7 populations in rows and the MEM scores for all 25 measured markers in columns. See \code{\link{PBMC}} dataset for more details.
}
\examples{
data(MEM_matrix)
}
\keyword{datasets}
