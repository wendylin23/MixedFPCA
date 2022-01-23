##########################################
######### Basis Generator ################
##########################################

library(fda)

### Fourier basis
f1.basis <- function(d, k)
{ 
  D = length(d)
  1/sqrt(D) * ((k%%2) * sin((k + k%%2) * pi * (d - 0.5)/D) + 
                 (1 - k%%2) * cos((k + k%%2) * pi * (d - 0.5)/D))
}

f2.basis <- function(d, k)
{
  D = length(d)
  if (k == 1) 
    rep(1/sqrt(2 * D), D)
  else
    1/sqrt(D) * ((k%%2) * cos((k + k%%2 + 4) * pi * (d - 0.5)/D) + 
                   (1 - k%%2) * sin((k + 4) * pi * (d - 0.5)/D))
}

### polynomial basis
p.basis <- function(d, K ,k){
  D = length(d)
  #as.function(legendre.polynomials(k, normalized = TRUE)[[k]])((2*(d - 0.5)/D - 1))/sqrt(D) * (sqrt(1/2))
  r <- chebyshev.t.recurrences(K, normalized=TRUE )
  as.function(orthogonal.polynomials( r )[[k+1]])((2*(d - 0.5)/D - 1))/sqrt(D) * (sqrt(1/2))
}

## error generator
err.generator <- function(nbasis,d,n,err.sd)
{
  D = length(d)
  basisobj = create.bspline.basis(c(1,D),nbasis)
  err.bs = smooth.basis(1:D, t(matrix(rnorm(n*D,0,err.sd), nrow = n)),basisobj)
  err.sm = t(eval.fd(1:D,err.bs$fd))
}

### plot basis function
plot_basis <- function(phi.mat,D, add = FALSE)
{
  if(add == FALSE)
  {
    plot(1, type="n", xlab="", ylab="", xlim=c(1,D), ylim=range(phi.mat))
    for (i in 1:dim(phi.mat)[1]) {
      lines(1:D, phi.mat[i,],col=i)
    } 
  }else{
    for (i in 1:dim(phi.mat)[1]) {
      lines(1:D, phi.mat[i,],col=i)
    }
  }
  #cormat <- apply(phi.mat, MARGIN=1, FUN=function(z) apply(phi.mat, MARGIN=1, FUN=function(y) cor(z, y)))
  #cormat
}

### simulation generator
MixedFPCA_simulation = function(params, levels, basis_list, times = NULL)
{
  M           <- params$M # - number of subjects 
  J           <- params$J
  S           <- params$S
  balanced    <- params$balanced          # balanced ("ba") or unbalanced ("ub") design           
  
  
}
