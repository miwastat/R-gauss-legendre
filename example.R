### Calsulates nodes and weights of Gauss-Legendre quadrature.
### The nodes and weights are separated by ",\n"
### so that they can be used in other programs (C and R).
###
### For 128-bit quadruple values (with 112-bit fraction part),
###   the precision is about log10(2^113) ~ 34 decimal places.
### So 36-digit constants are printed.
###

library(Rmpfr)
source("gauss-legendre.R")

bits <- 160   # log10(2^160) ~ 48 decimal places
digits <- 36  # 36 decimal digits for printing
eps <- mpfr(10, bits)^(-42) # precision for Newton procedure

n <- 20       # number of nodes for Gauss-Legendre quadrature
y <- gauss.legendre(n, bits=bits, eps=eps)

## Print nodes and weights for copying & pasting into other programs.
cat(formatMpfr(y$nodes, digits=digits), sep=",\n")
cat(formatMpfr(y$weights, digits=digits), sep=",\n")

## Scientific expressions
cat(formatMpfr(y$nodes, digits=digits, scientific=TRUE), sep=",\n")
cat(formatMpfr(y$weights, digits=digits, scientific=TRUE), sep=",\n")
