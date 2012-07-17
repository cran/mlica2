\name{mlica}
\alias{mlica}

\title{Maximum likelihood implementation of independent component analysis}

\description{This function performs ICA using a maximum likelihood framework
     and takes as arguments parameters to control the number of
     algorithm runs and convergence criteria.}
\usage{
mlica(prNCP, nruns = 10, tol = 1e-04, maxit = 300, fail.th = 5, learn.mu = 1)
}

\arguments{
  \item{prNCP}{The output object from 'proposeNCP'.}
  \item{nruns}{The number of converged algorithm runs sought (function
          returns the best solution according to the log-likelihood
          value).}
  \item{tol}{Tolerance level for establishing convergence of run.}
  \item{maxit}{Maximum number of iterations to allow per run.}
  \item{fail.th}{A threshold on the number of consecutive runs that fail to
          converge.}
  \item{learn.mu}{Learning parameter for fixed point algorithm (note that this
          need not be changed since it has already been optimised).}
}

\value{
A list with following components: 

       A: Estimate of the mixing matrix.

       B: Estimate of the inverse mixing matrix.

       S: Estimate of the source matrix.

       X: Normalised data matrix.

     ncp: Number of independent components.

      NC: Binary number specifying whether best run converged or
          not.(=1 indicates convergence,=0 indicates no convergence).

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

	Chiappetta P., Roubaud MC. and Torresani B.: Blind source separation
     and the analysis of microarray data, J. Comput. Biol. 2004;
     11(6):1090-109.
}
\author{Andrew Teschendorff a.teschendorff@ucl.ac.uk}

\examples{
\dontrun{
data(simMAdata);
dataX <- simMAdata[[1]];
prPCA <- PriorNormPCA(dataX);
prNCP <- proposeNCP(prPCA,0.1);
a.best.l <- list();
for( i in 1:5){
 a.best.l[[i]] <- mlica(prNCP,nruns=5);
}
checkICA <- CheckStability(a.best.l,0.7);
sourceS <- simMAdata[[3]];
print(cor(a.best.l[[1]]$S,sourceS));
sModes <- SortModes(a.best.l[[1]],c.val=0.5);
}

## The function is currently defined as
function (prNCP, nruns = 10, tol = 1e-04, maxit = 300, fail.th = 5, 
    learn.mu = 1) 
{
    print("Entering mlica")
    print("Performing preliminary run")
    a <- mlicaMAIN(prNCP, tol = 1e-04, maxit = 10, mu = learn.mu)
    ncp <- dim(a$S)[2]
    max.logL <- a$LL
    a.best <- a
    print("Finished preliminary run")
    print("Starting runs")
    run.n <- 0
    fail.count <- 0
    v.logL <- vector()
    v.NC <- vector()
    while (run.n < nruns) {
        a <- mlicaMAIN(prNCP, tol = 1e-04, maxit = maxit, mu = learn.mu)
        v.logL <- c(v.logL, a$LL)
        v.NC <- c(v.NC, a$NC)
        if (a$NC == 0) {
            fail.count <- 0
            run.n <- run.n + 1
            if (a$LL > max.logL) {
                a.best <- a
            }
        }
        else {
            fail.count <- fail.count + 1
        }
        if (fail.count >= fail.th) {
            print("Stopping: Five consecutive runs failed to converge!")
            print("Consider either increasing the threshold for pca eigenvalues to perform ICA on a smaller subspace or increasing maxit")
            stop
        }
    }
    print("End of runs")
    return(a.best)
  }
}
