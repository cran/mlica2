\name{mlicaMAIN}
\alias{mlicaMAIN}

\title{Main engine function that implements the fixed point algorithm for maximum likelihood of ICA modes}
\description{See references for detailed description.}
\usage{
mlicaMAIN(prNCP, tol = 1e-04, maxit = 300, mu = 1)
}

\arguments{
  \item{prNCP}{The output object of 'proposeNCP'.}
  \item{tol}{Tolerance level for convergence.}
  \item{maxit}{Maximum number of iterations to allow for convergence.}
  \item{mu}{Learning paramter for fixed point algorithm. This has already
          been optimised.}
}

\value{A list with following components: 

       A: Estimate of the mixing matrix.

       B: Estimate of the inverse mixing matrix.

       S: Estimate of the source matrix.

       X: Normalised data matrix.

     ncp: Number of independent components.

      NC: Binary number specifying whether best run converged,0, or
          not,1.

      LL: Log likelihood value of best run.
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
function (prNCP, tol = 1e-04, maxit = 300, mu = 1) 
{
    print("Entered MLica")
    X <- prNCP$X
    x <- prNCP$x
    pEx <- prNCP$pEx
    pCorr <- prNCP$pCorr
    ntp <- dim(X)[1]
    ndim <- dim(X)[2]
    ncp <- ncol(x)
    Sest <- matrix(nrow = ntp, ncol = ncp)
    B.old <- matrix(runif(ncp * ncp, 0, 1), nrow = ncp, ncol = ncp)
    B.o <- B.old
    icount <- 0
    not.conv <- c(1, 2)
    y <- matrix(nrow = ntp, ncol = ncp)
    tmp <- matrix(nrow = ncp, ncol = ncp)
    beta <- vector(length = ncp)
    alpha <- vector(length = ncp)
    while ((length(not.conv) > 0) && (icount < maxit)) {
        print(c("Entering iteration loop ", icount))
        Cy <- B.old \%*\% t(B.old)
        svds <- eigen(Cy, symmetric = TRUE)
        D <- diag(svds$values)
        E <- svds$vectors
        Dinv <- solve(D)
        V <- E \%*\% sqrt(Dinv) \%*\% t(E)
        B.old <- V \%*\% B.old
        for (g in 1:ntp) {
            y[g, ] <- B.old \%*\% x[g, ]
        }
        for (c in 1:ncp) {
            beta[c] <- 2 * sum(y[, c] * tanh(y[, c]))/ntp
            alpha[c] <- -1/(beta[c] - 2 + 2 * sum(tanh(y[, c]) * 
                tanh(y[, c]))/ntp)
            for (c2 in 1:ncp) {
                tmp[c, c2] <- -2 * sum(tanh(y[, c]) * y[, c2])/ntp
            }
        }
        print("Checkpt1")
        tmp <- diag(beta) + tmp
        B <- B.old + mu * diag(alpha) \%*\% tmp \%*\% B.old
        Dev <- abs(B - B.o)
        AvDev <- sum(Dev)/(ncp * ncp)
        print(c("AvDev=", AvDev))
        not.conv <- vector()
        not.conv <- as.vector(Dev[Dev > tol])
        B.old <- B
        B.o <- B
        icount <- icount + 1
        for (g in 1:ntp) {
            Sest[g, ] <- B \%*\% x[g, ]
        }
        logL <- -2 * sum(log(cosh(Sest))) + ntp * log(abs(det(B)))
        print("iterated logL")
        print(logL)
    }
    Cy <- B \%*\% t(B)
    svds <- eigen(Cy, symmetric = TRUE)
    D <- diag(svds$values)
    E <- svds$vectors
    Dinv <- solve(D)
    V <- E \%*\% sqrt(Dinv) \%*\% t(E)
    B <- V \%*\% B
    for (g in 1:ntp) {
        Sest[g, ] <- B \%*\% x[g, ]
    }
    Aest <- t(pEx \%*\% sqrt(pCorr) \%*\% t(B))
    if (length(not.conv) > 0) {
        NotConv <- 1
    }
    else {
        NotConv <- 0
    }
    logL <- -2 * sum(log(cosh(Sest))) + ntp * log(abs(det(B)))
    return(list(A = Aest, B = B, S = Sest, X = X, ncp = dim(Sest)[2], 
        NC = NotConv, LL = logL))
  }
}
