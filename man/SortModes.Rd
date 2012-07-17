\name{SortModes}
\alias{SortModes}

\title{Sorting of ICA Modes}
\description{Sorts inferred ICA modes using two criteria: Relative data power
     or the Liebermeister criterion, which is based on a measure that
     is a weighted linear combination of non-gaussianity and data
     variance measures.}
\usage{
SortModes(a.best,c.val = 0.25)
}

\arguments{
  \item{a.best}{The output object of 'mlica'.}
  \item{c.val}{A parameter to control the relative weight of the two
          measures when using the Liebermeister criterion. Should be
          between 0 (pure data variance measure) and 1 (pure
          non-gaussianity).}
}

\value{A list with the components:

	a.best: The output of 'mlica'.

	rdp: The relative data power values obtained for each independent
          component.

	lbm: The Liebermeister contrast value for each component.
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

#This function is currently defined as
function(a.best,c.val=0.25){
  
	ncp <- ncol(a.best$S);
	Ng <- nrow(a.best$S);

	#SORTING CRITERION

	# A) Computation of relative data power. Store values in vector of size H=ncp. Could use different criterion here.
	# need squared entries
	 Ssq <- a.best$S * a.best$S ; 
	 Asq <- a.best$A * a.best$A ;
	 Xsq <- a.best$X * a.best$X ;
	 rdp <- rep(0, times=ncp);
	 for ( k in 1:ncp ) {
	  rdp[k]<- sum(Ssq[,k])*sum(Asq[k,])/sum(Xsq) ;
	 }
	 rdp.s <- sort(rdp, na.last=NA,decreasing=TRUE, index.return=TRUE);

	# B) sorting with mixture of contrast and data variance (Liebermeister)
	 JG <- rep(0, times=ncp);
	 JA <- rep(0, times=ncp);
	 # generate values from std. normal distribution
	 nu <- rnorm(10000,0,1);
	 G0 <- mean(log(cosh(nu)));
	 
	 for ( k in 1:ncp ){
	   # compute contrast for mode using logcosh
	   G1 <- mean(log(cosh(a.best$S[,k])));
	   JG[k] <- abs(G1-G0);
	   JA[k] <- sum(Asq[k,]);
	 }   
	   sumJG <- sum(JG) ; sumJA <- sum(JA) ;
	   J <- JG*(c.val/sumJG)+ JA*(1-c.val)/sumJA ;
	   J.s <- sort(J, na.last=NA,decreasing=TRUE, index.return=TRUE);

	   return(list(a.best=a.best,rdp=rdp.s,lbm=J.s));

	} 
}

\keyword{ correlation }
\keyword{ ICA }
