\name{CheckStability}
\alias{CheckStability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test stability of inferred ICA modes}
\description{Performs a correlation test to see which of the inferred ICA modes
     are reproducible across multiple runs using different random
     initialisations. Returns a set of consensus ICA modes and
     stability scores for each following the algorithm of
     Chiappetta,...et.al (2004)}
\usage{
CheckStability(a.best.l, corr.th = 0.7)
}

\arguments{
  \item{a.best.l}{List of 'a.best' objects from 'mlica' runs.}
  \item{corr.th}{Correlation threshold to use to decide whether a mode is
          reproducible.}
}

\value{A list with the following components 

   consS: Consensus source matrix with columns labeling the consensus
          ICA modes. Has same number of rows as 'a.best$S'.

   consA: Consensus mixing matrix with rows labeling the consensus ICA
          modes.

   stabM: Vector of same length as 'consM' giving the stability
          measures of each consensus ICA mode. Stability or
          reproducibility measures are given as fractions, that is, the
          number of times the ICA mode correlates with one of the other
          runs at threshold level 'corr.th' divided by the number of
          runs (length of 'a.best.l').}
\references{
	Hyvaerinen A., Karhunen J., and Oja E.: Independent Component
     Analysis, John Wiley and Sons, New York, (2001).

	Kreil D. and MacKay D. (2003): Reproducibility Assessment of
     Independent Component Analysis of Expression Ratios from DNA
     microarrays, Comparative and Functional Genomics *4* (3),300-317.

	Liebermeister W. (2002): Linear Modes of gene expression determined
     by independent component analysis, Bioinformatics *18*, no.1,
     51-60.

	Chiappetta P., Roubaud MC. and Torresani B.: Blind source separation
     and the analysis of microarray data, J. Comput. Biol. 2004;
     11(6):1090-109.
}
\author{Andrew Teschendorff a.teschendorff@ucl.ac.uk}

\examples{

## The function is currently defined as
function (a.best.l, corr.th = 0.7) 
{
    nIruns <- length(a.best.l)
    if (nIruns < 2) {
        print("Argument must be list of at least two a.best objects")
        stop
    }
    for (i in 2:nIruns) {
        if (a.best.l[[1]]$ncp != a.best.l[[i]]$ncp) {
            print("Stopping: the objects in list must have same number of ICA modes")
            stop
        }
    }
    ncp <- a.best.l[[1]]$ncp
    X <- a.best.l[[1]]$X
    nConn.v <- vector(length = nIruns * ncp)
    Adj.a <- array(rep(0, nIruns * ncp * nIruns * ncp), dim = c(nIruns, 
        ncp, nIruns, ncp))
    v <- 1
    for (j1 in 1:nIruns) {
        selrun <- setdiff(1:nIruns, j1)
        for (j2 in selrun) {
            for (i1 in 1:ncp) {
                cor.v <- as.vector(abs(cor(a.best.l[[j1]]$S[, 
                  i1], a.best.l[[j2]]$S)))
                i2 <- which.max(cor.v)
                Adj.a[j1, i1, j2, i2] <- 1
            }
        }
        for (i in 1:ncp) {
            nConn <- 0
            for (j2 in selrun) {
                if (sum(Adj.a[j1, i, j2, ]) > 0) {
                  nConn <- nConn + 1
                }
            }
            nConn.v[v] <- nConn
            v <- v + 1
        }
    }
    selmodes.idx <- 1:(nIruns * ncp)
    consS.lv <- list()
    StabScore.l <- list()
    c <- 1
    while (c <= ncp) {
        max.idx <- which.max(nConn.v)
        j1.max <- as.integer(max.idx/(ncp + 1)) + 1
        i1.max <- max.idx - (j1.max - 1) * ncp
        selrun <- setdiff(1:nIruns, j1.max)
        consS <- a.best.l[[j1.max]]$S[, i1.max]
        corr.idx <- vector()
        for (j2 in selrun) {
            i2 <- which(Adj.a[j1.max, i1.max, j2, ] == 1)
            if (length(i2) > 0) {
                corr.idx <- c(corr.idx, (j2 - 1) * ncp + i2)
                tmp <- cor(a.best.l[[j1.max]]$S[, i1.max], a.best.l[[j2]]$S[, 
                  i2])
                consS <- consS + sign(tmp) * a.best.l[[j2]]$S[, 
                  i2]
            }
        }
        selmodes.idx <- setdiff(selmodes.idx, c(max.idx, corr.idx))
        consS.lv[[c]] <- consS/(nConn.v[max.idx] + 1)
        StabScore.l[[c]] <- (nConn.v[max.idx] + 1)/nIruns
        nConn.v[c(max.idx, corr.idx)] <- -1
        c <- c + 1
    }
    consS.m <- matrix(nrow = nrow(a.best.l[[1]]$S), ncol = ncp)
    StabScore.v <- vector()
    for (c in 1:ncol(consS.m)) {
        consS.m[, c] <- consS.lv[[c]]/sqrt(var(consS.lv[[c]]))
        StabScore.v[c] <- StabScore.l[[c]]
    }
    STS <- t(consS.m) \%*\% consS.m
    STSinv <- solve(STS)
    STX <- t(consS.m) \%*\% X
    consA.m <- STSinv \%*\% STX
    return(list(stabM = StabScore.v, consS = consS.m, consA = consA.m))
  }
}

\keyword{ correlation }
\keyword{ ICA }
