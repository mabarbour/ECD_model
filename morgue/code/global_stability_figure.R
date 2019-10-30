## load required data ----
C1.ESS.df <- read.csv('data/evol.symmetry.df.C1.3sp.csv')
sym.4sp.df <- read.csv('data/evol.symmetry.df.4sp.csv')

## load required functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

#### Setup the initial state variables and global parameters for all of the models ----

## set state variables for all models
# These are different from the equilbrium values to illustrate the dynamic trajectories back to equilibrium.
R1 <- 2
R2 <- 2
C1 <- 1
C2 <- 1
i.state.2C_2R <- c(R1 = R1, R2 = R2, C1 = C1, C2 = C2) 

## set duration of simulations
Time <- 1000

## set up general code for the plot. Maintaining 'w' as the key parameter to manipulate
#aii <- C1.ESS.df$a11[1]
#aij <- C1.ESS.df$a12[1]

#aii_seq <- seq(aii, aii*2, by = 0.01)[1:75]
#aij_seq <- seq(aij, 0, by = -0.01)[1:75]

#w11 <- w22 <- C1.ESS.df$w11[1]
#a_alt <- a_start
#aECD <- w11*aii_seq + (1-w11)*aij_seq

# general parameters. Taken from Figure 3 in McCann et al. 2005. 
#r1 <- r2 <- 1
#K1 <- K2 <- 3.25
#e11 <- e12 <- e21 <- e22 <- 0.8
#h11 <- h12 <- h21 <- h22 <- 0.4
#m1 <- m2 <- 1
#w11 <- w22 <- C1.ESS.df$w11[1]

#aii <- C1.ESS.df$a11[1]
#aij <- C1.ESS.df$a12[1]

#aii_seq <- seq(aii, aii*2, by = 0.01)[1:75]
#aij_seq <- seq(aij, 0, by = -0.01)[1:75]

p.set <- 1

p1 <- sym.4sp.df[p.set, 2:21]

#param.1 <- c(r1 = r1, r2 = r2, K1 = K1, K2 = K2,
#             e11 = e11, e12 = e12, e21 = e21, e22 = e22,
#             h11 = h11, h12 = h12, h21 = h21, h22 = h22,
#             a11 = aii_seq[1], a12 = aij_seq[1], a21 = aij_seq[1], a22 = aii_seq[1],
#             m1 = m1, m2 = m2, w11 = w11, w22 = w22)

# Run the experiment.
init <- ode(i.state.2C_2R, 1:Time, ECD_model, p1)

## effective attack rate. A key component in this model
aEFF.1 <- p1[1,"a11"]*p1[1,"w11"] + p1[1,"a12"]*(1 - p1[1,"w11"])
aEFF.2 <- p1[1,"a22"]*p1[1,"w22"] + p1[1,"a21"]*(1 - p1[1,"w22"])

## Resource 1 isocline
Rx <- seq(0.1, p1[1,"K1"], 0.1) # manipulating different Resource densities to solve R isocline.

# set R1 = 0 and solved algebraically. Note R1 = R2.
Riso <- expression(
  (p1[1,"r1"]*(p1[1,"K1"] - Rx)*(1 + p1[1,"h11"]*Rx*aEFF.1)) / 
  (p1[1,"K1"]*aEFF.1) +
    (p1[1,"r2"]*(p1[1,"K2"] - Rx)*(1 + p1[1,"h22"]*Rx*aEFF.2)) / 
    (p1[1,"K2"]*aEFF.2)) 

RisoStable <- eval(Riso)

## Consumer isocline

# set C1 = 0, and solved algebraically. Note C1 = C2.
Ciso <- expression(
  p1[1,"m1"] / 
  ((p1[1,"e11"] - p1[1,"h11"] * p1[1,"m1"]) * aEFF.1) +
    p1[1,"m2"] / 
    ((p1[1,"e22"] - p1[1,"h22"] * p1[1,"m2"]) * aEFF.2)) 

CisoStable <- eval(Ciso)

## Plot stability around consumer and resource isoclines
plot(Rx,RisoStable, type = "l", ylim = c(0,3), xlim = c(0,5), ylab = "Consumer density", xlab = "Resource density")
abline(v=CisoStable, lty = 2, col =2) 
legend("topleft", c("R-isocline","C-isocline"), lty=1:2, bty="n", cex=0.8, col=1:2)

points(I(i.state.2C_2R[1] + i.state.2C_2R[2]),
       I(i.state.2C_2R[3] + i.state.2C_2R[4])) # starting point of experiment
arrows(x0 = I(init[-Time,2]+init[-Time,3]), y0 = I(init[-Time,4]+init[-Time,5]), 
       x1 = I(init[-1,2]+init[-1,3]), y1 = I(init[-1,4]+init[-1,5]), length=0.1, lty=1) # trace stability across different time steps. Don't know why there are so many warnings

test <- filter(sym.4sp.df, mut.suc != 0)
dim(test)
plot(-1*max.Re.eigen ~ c(1:dim(test)[1]), test)
plot(I(a11*w11 + a12*(1-w11)) ~ c(1:dim(test)[1]), test)

## Second set of parameters
p1 <- sym.4sp.df[1, 2:21]

param.2 <- c(r1 = r1, r2 = r2, K1 = K1, K2 = K2,
             e11 = e11, e12 = e12, e21 = e21, e22 = e22,
             h11 = h11, h12 = h12, h21 = h21, h22 = h22,
             a11 = aii_seq[75], a12 = aij_seq[75], a21 = aij_seq[75], a22 = aii_seq[75],
             m1 = m1, m2 = m2, w11 = w11, w22 = w22)

# Run the experiment.
init.2 <- ode(i.state.2C_2R, 1:Time, ECD_model, param.2)

## effective attack rate. A key component in this model
aEFF <- aii_seq[75]*w11 + aij_seq[75]*(1 - w11)

## Resource isocline
Rx <- seq(0.1,K1,0.1) # manipulating different Resource densities to solve R isocline.
Riso <- expression((r1*(K1-Rx)*(1 + h11*Rx*aEFF))/(K1*aEFF)) # set R1 = 0 and solved algebraically. Note R1 = R2.
RisoStable <- eval(Riso)

## Consumer isocline
Ciso <- expression(m1 / ((e11 - h11*m1)*aEFF)) # set C1 = 0, and solved algebraically. Note C1 = C2.
CisoStable <- eval(Ciso)

## Plot stability around consumer and resource isoclines
plot(Rx,RisoStable, type = "l", ylim = c(0,2), ylab = "Consumer density", xlab = "Resource density")
abline(v=CisoStable, lty = 2, col =2) 
legend("topleft", c("R-isocline","C-isocline"), lty=1:2, bty="n", cex=0.8, col=1:2)
points(i.state.2C_2R[1],i.state.2C_2R[3]) # starting point of experiment
arrows(x0 = init.2[-Time,2], y0 = init.2[-Time,4], x1 = init.2[-1,2], y1 = init.2[-1,4], length=0.1, lty=1) # trace stability across different time steps. Don't know why there are so many warnings



