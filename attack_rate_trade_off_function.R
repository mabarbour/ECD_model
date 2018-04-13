# from Sargent and Otto 2006
# (a11/A)^n + (a12/A)^n = 1

## using Mathematica
a_ij_tradeoff <- function(A, a_ii, n) {
  A * (1 - (a_ii/A)^n)^(1/n)
}
