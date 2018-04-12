# from Sargent and Otto 2006
# (a11/A)^n + (a12/A)^n = 1

## using Mathematica
a_ij_tradeoff <- function(A, a_ii, n) {
  A * (1 - (a_ii/A)^n)^(1/n)
}

a11_seq <- seq(0,4,0.01)

a_ij_tradeoff(A=4, a_ii=a11_seq, n=1.5)

## CONFIRM TRADE-OFF
plot(a_ij_tradeoff(A=4, a_ii=a11_seq, n=1) ~ a11_seq, type="l")
lines(a_ij_tradeoff(A=4, a_ii=a11_seq, n=1.1) ~ a11_seq, type="l", col="blue")
lines(a_ij_tradeoff(A=4, a_ii=a11_seq, n=0.9) ~ a11_seq, type="l", col="red")
