\name{proposeNCP}
\alias{proposeNCP}

\title{Number of independent components proposal function}
\description{This function takes the output of 'PriorNormPCA' and returns for a
     given threshold the number of components to be inferred for
     subsequent ICA.}
\usage{
proposeNCP(prPCA, thresh = 0.1)
}

\arguments{
  \item{prPCA}{The output object from 'PriorNormPCA'.}
  \item{thresh}{Threshold on eigenvalues.}
}

\value{A list with following components: 

       X: Normalised data matrix.

       x: Normalised data matrix projected onto selected subspace.

     pEx: Selected eigenvectors defining subspace for projection.

   pCorr: Projected correlation matrix.

     ncp: The dimension of the selected subspace(=number of independent
          components to be inferred with subsequent ICA).}
\references{
	Hyvaerinen A., Karhunen J., and Oja E.: Independent Component
     Analysis, John Wiley and Sons, New York, (2001).

	Kreil D. and MacKay D. (2003): Reproducibility Assessment of
     Independent Component Analysis of Expression Ratios from DNA
     microarrays, Comparative and Functional Genomics *4* (3),300-317.

	Liebermeister W. (2002): Linear Modes of gene expression determined
     by independent component analysis, Bioinformatics *18*, no.1,
     51-60.
}
\author{Andrew Teschendorff a.teschendorff@ucl.ac.uk}

\examples{

## The function is currently defined as
function (prPCA, thresh = 0.1) 
{
    X <- prPCA$X
    eigenvals.v <- diag(prPCA$Dx)
    Ex <- prPCA$Ex
    ntp <- nrow(X)
    ndim <- ncol(X)
    print("About to find ncp")
    p.cpts <- eigenvals.v[eigenvals.v > thresh]
    ncp <- length(p.cpts)
    pCorr <- diag(eigenvals.v[1:ncp])
    pEx <- Ex[, 1:ncp]
    x <- matrix(nrow = ntp, ncol = ncp)
    for (g in 1:ntp) {
        for (c in 1:ncp) {
            x[g, c] <- sum(X[g, ] * Ex[, c])/sqrt(diag(pCorr)[c])
        }
    }
    return(list(X = X, x = x, pEx = pEx, pCorr = pCorr, ncp = ncp))
  }
}
