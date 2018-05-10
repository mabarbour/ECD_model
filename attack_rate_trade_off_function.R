# from Sargent and Otto 2006
# (a11/A)^n + (a12/A)^n = 1

## using Mathematica
a_ij_tradeoff <- function(A, a_ii, n) {
  A * (1 - (a_ii/A)^n)^(1/n)
}
a_ij_tradeoff(A=2, a_ii=1.5, n=1)

# from Abrams 1986
# ki*C1i^2 + (1-ki)*C2i^2 = ki*(1-ki)
a_ij_Abrams <- function(trait, k, n) {
  ((k * (1 - k) - k * trait^n)/(1 - k))^(1/n)
}
df <- data.frame(B11 = 1)
a_ij_Abrams(trait = 1, k=0.5, n=1)

a_ii_seq <- seq(0,2,0.01)

plot(a_ij_Abrams(ki=0.7, a_ii=a_ii_seq, n=1) ~ a_ii_seq, type='l')
