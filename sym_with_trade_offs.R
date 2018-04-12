# from Sargent and Otto 2006
# (a11/A)^n + (a12/A)^n = 1
# (a12/A)^n = 1 - (a11/A)^n    
# (a12/A) = 1^(1/n) - (a11/A)
# a12 = A*1^(1/n) - a11

library(rootSolve)

# I can explore MacArthur model dynamics by setting h=0 and w=0
# Explore linear and concave trade-off scenarios. No evidence for convex one in the stickleback system.

## need to get trade-off calculation correct for 

# note that a11 must be greater than half of A for the function to work
a_trade_off <- function(A, a11, a12, n) (a11/A)^n + (a12/A)^n - 1

a12_root <- uniroot(f = a_trade_off, interval = c(0,2), tol = 0.001, A=4, a11=1, n=-0.5)$root

a11_seq <- seq(2,4,0.1)
a12_roots <- c()
for(i in 1:length(a11_seq)){
  a12_roots[i] <- uniroot(f = a_trade_off, interval = c(0,2), tol = 0.001, A=4, a11=a11_seq[i], n=-0.5)$root
}

uniroot.all(a_trade_off(A = 5, a11=a11, n=1), c(0,5))


a12_lin <- a_trade_off(A = 5, a11 = a11, n=1)
plot(a12_lin, a11)

a12_cc <- a_trade_off(A = 5, a11 = a11, n=-1)
plot(a12_cc, a11)

plot(a12_lin, a11, type='l')
lines(a12_cc, a11)
