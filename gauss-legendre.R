### gauss.legendre(n, bits, eps)
###   returns nodes and weights with arbitrary precision
###   for the Gauss-Legendre quadrature.
###
###  Arguments
###   n:    number of nodes
###   bits: number of bits for mpfr values in Rmpfr package
###   eps:  precision for the root of P(x, n)=0
###
### Output
###   gauss.legendre()$nodes:   nodes
###   gauss.legendre()$weights: weights
###   gauss.legendre()$maxitr:  max number of iterations
###
### Note
###   1) Rmpfr should be libraried before this file is sourced.
###      See 'example.R'.
###   2) Rmpfr's license is GPLv2 (or later).
###      Rmpfr interfaces to 'MPFR' with 'LGPL' license.
###
### History
###   2018-09-23: First coded.
###
### License
###   GPLv3 (Free and No warranty)
###   https://www.gnu.org/licenses/
###
### Coded by Tetsuhisa Miwa
###


### n-th order Legendre polynomial
P <- function(x, n, bits=160) {
  ## returns the n-th order Legendre polynomial value.
  ## x should be an mpfr number (-1.0 < x < 1.0 for P(x, n)=0).
  ## This function returns 1st and 2nd derivarives and P(x, n-1).

  p0 <- mpfr(1, bits)
  p1 <- x
  for(r in 2:n) {
    p  <- ((2*r-1)*x*p1 - (r-1)*p0)/r
    p0 <- p1
    p1 <- p
  }

  ## f=P(x, n), fd=P'(x, n), fdd=P''(x, n), f1=P(x, n-1)
  ## f1=P(x, n-1) is used for calculating weights.
  fd <- n*(p0-x*p1)/(1-x*x)
  fdd <- (2*x*fd - n*(n+1)*p1)/(1-x*x)
  return(list(f=p1, fd=fd, fdd=fdd, f1=p0))
}

### Newton method to solve x such that P(x, n)=0
newton <- function(x, n, bits=160, eps=mpfr(10, bits)^(-42)) {
  ## The 2nd-order Newton method is used.
  ## x:   initial value
  ## eps: precision both for x and P(x, n)

  delta <- mpfr(0, bits)
  for(i in 1:100){
    p.n <- P(x, n=n, bits=bits)
    f <- p.n$f
    if((abs(delta) < eps) && (abs(f) < eps))
      break
    fd <- p.n$fd            # 1st derivative
    fdd <- p.n$fdd          # 2nd derivative
    dscr <- fd*fd - 2*f*fdd # discriminant of 2nd-order polynomial
    if(dscr > 0) {
      if(fd > 0) {
        delta <- 2*f / (-fd-sqrt(dscr))
      } else {
        delta <- 2*f / (-fd+sqrt(dscr))
      }
    } else {
      if(abs(fdd) < eps || abs(fd) < eps) {
        delta <- eps  # Perturb if we cannot move.
      } else {
        delta <- -fd/fdd
      }
    }
    x <- x + delta
  }
  ## The procedure converged if itr < 100.
  return(list(x=x, f1=P(x, n, bits)$f1, itr=i))
}

### Nodes and weights for Gauss-Legendre quadrature
gauss.legendre <- function(n, bits=160, eps=mpfr(10, bits)^(-42)) {
  np <- (n+1)%/%2   # np=n/2 for even n, and np=(n+1)/2 for odd n
  pi.mpfr <- Const("pi", bits)
  maxitr <- 0       # max number of iterations
  nodes <- mpfrArray(NA, precBits=bits, dim=np)
  weights <- mpfrArray(NA, precBits=bits, dim=np)
  
  for(k in 1:np) {
    x <- sin(pi.mpfr*(n+1-2*k)/(2*n+1))  # initial approximation
    zero <- newton(x, n=n, bits=bits, eps=eps)
    nodes[k] <- zero$x
    weights[k] <- 2*(1-(zero$x)*(zero$x))/(n*zero$f1)^2
    ## Comment out the next line if you need not interim messages.
    print(sprintf("Calculating k=%i, iterations=%i", k, zero$itr))
  }
  return(list(nodes=nodes, weights=weights, maxitr=maxitr))
}
