\name{PriorNormPCA}
\alias{PriorNormPCA}

\title{Prior PCA analysis for threshold setting and noise removal}
\description{This function performs a simple PCA analysis to aid in threshold
     setting and noise removal.}
\usage{
PriorNormPCA(X)
}

\arguments{
  \item{X}{Data Matrix (need not be normalised). Subsequent ICA seeks
          independent modes as independent distributions with values
          "down the rows".}
}
\details{This function performs a simple PCA analysis and is used prior to
     application of the main ICA algorithm. The objective of the prior
     PCA is to help determine the dimensionality of a subspace on which
     the further ICA converges. The convention used here is that the
     rows of 'X' label the space over which independent components are
     sought. For a typical microarray application in which ICA is being
     used as a generative model for gene expression, rows should label
     genes and columns should label samples. If, however, ICA is to be
     used as an unsupervised projection pursuit algorithm, rows should
     label samples and columns genes. For the latter application, the
     number of genes should be less than the number of samples.}
\value{A list with following components: 

       X: Normalised data matrix with the mean of each column set to
          zero.

      Dx: Eigenvalues in a diagonal matrix.

      Ex: Eigenvectors
}
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
function (X) 
{
    ndim <- ncol(X)
    ntp <- nrow(X)
    for (s in 1:ndim) {
        X[, s] <- X[, s] - mean(X[, s])
    }
    print("Performing SVD")
    svd.o <- svd(X, LINPACK = TRUE)
    Dx <- diag(svd.o$d * svd.o$d)/ntp
    Ex <- svd.o$v
    barplot(Dx, main = "Singular values")
    return(list(X = X, Dx = Dx, Ex = Ex))
  }
}