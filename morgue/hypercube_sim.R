######################################################
## Description: Identify feasible and stable values for parameters and state variables
## Code author(s): Matt Barbour
## Email: barbour@zoology.ubc.ca
######################################################

## load required libraries ----
library(deSolve)
library(rootSolve)
library(plyr)
library(dplyr)
library(assertr)
library(pse) # for latin hypercube

## load required functions ----
source('Google Drive/ms_ecd/numerical_fxns.R')
source('Google Drive/ms_ecd/local_stability.R')

## setup latin hypercube sampling for parameter sets ----
n.p <- 30 # number of parameter sets
factors <- c("R1_0","R2_0","C1_0","C2_0",
             "r1","r2","K1","K2",
             "e11","e12","e21","e22",
             "h11","h12","h21","h22",
             "a11","a12","a21","a22",
             "m1","m2","w11","w22")
q <- replicate(24, "qunif")
q.arg <- list(list(min= 0.5, max = 2),list(min = 0.5, max = 2),list(min = 0.5, max = 2),list(min = 0.5, max = 2),
           list(min = 0.5, max = 2),list(min = 0.5, max = 2),list(min = 2, max = 5),list(min = 2, max = 5),
           list(min = 0.5, max = 0.9),list(min = 0.5, max = 0.9),list(min = 0.5, max = 0.9),list(min = 0.5, max = 0.9),
           list(min = 0.1, max = 0.6),list(min = 0.1, max = 0.6),list(min = 0.1, max = 0.6),list(min = 0.1, max = 0.6),
           list(min = 0.1, max = 4),list(min = 0.1, max = 4),list(min = 0.1, max = 4),list(min = 0.1, max = 4),
           list(min = 0.1, max = 1),list(min = 0.1, max = 1),list(min = 0.05, max = 0.95),list(min = 0.05, max = 0.95))

LHS.4 <- LHS(model=NULL, factors=factors, N=n.p, q=q, q.arg=q.arg, nboot=0) # use nboot=0 if do not need partial correlations

LHS.4.ps <- get.data(LHS.4) 

LHS.4.ps[1,5:24] <- c(r1 = 1, r2 = 1, K1 = 3, K2 = 3, e11 = 0.6, e12 = 0.6, e21 = 0.6, e22 = 0.6, h11 = 0.4, h12 = 0.4, h21 = 0.4, h22 = 0.4, a11 = 5, a12 = 0.5, a21 = 0.5, a22 = 5, m1 = 1, m2 = 1, w11 = 0.6, w22 = 0.6)

runsteady(parms = LHS.4.ps[3,5:24], y = c(R1 = LHS.4.ps[3,"R1_0"], R2 = LHS.4.ps[3,"R2_0"], C1 = LHS.4.ps[3,"C1_0"], C2 = LHS.4.ps[3,"C2_0"]), times = c(0,Inf), func = ECD_model)

tout <- failwith(default = c(R1 = NA, R2 = 2), runsteady)
test <- tout(parms = LHS.4.ps[1,5:24], y = c(R1 = LHS.4.ps[1,"R1_0"], R2 = LHS.4.ps[1,"R2_0"], C1 = LHS.4.ps[1,"C1_0"], C2 = LHS.4.ps[1,"C2_0"]), times = c(0,Inf), func = ECD_model)
class(test)
out$y
attr(out,"steady")
c(test$y, steady = attr(test,"steady"))

## generate for-loop for 4-species community ----

# Run simulation and create list with final state variables and whether they are at a steady state.

sim.length <- 1000 # duration of simulations on each parameter set
s.runsteady <- failwith(default = c(R1 = NA, R2 = NA, C1 = NA, C2 = NA, steady = FALSE), runsteady) # modify runsteady to still generate output with fatal error

list.LHS.4.ps <- list()
for(i in 1:dim(LHS.4.ps)[1]){
  out <- s.runsteady(parms = LHS.4.ps[i,5:24], 
                     y = c(R1 = LHS.4.ps[i,"R1_0"],
                           R2 = LHS.4.ps[i,"R2_0"],
                           C1 = LHS.4.ps[i,"C1_0"],
                           C2 = LHS.4.ps[i,"C2_0"]),
                     times = c(0,sim.length), func = ECD_model)
  
  
  
  out.df <- data.frame(R1 = out$y["R1"],
                       R2 = out$y["R2"],
                       C1 = out$y["C1"],
                       C2 = out$y["C2"],
                       steady = attr(out, "steady"))
  browser()
}

## Turn list into a data frame ----
df.4sp <- ldply(list.LHS.4.ps) %>% tbl_df()
as.matrix(df.4sp) # for viewing
#write.csv(x = df.4sp, file = "Google Drive/ms_ecd/df.4sp.csv")

## generate for-loop for 3-species community ----
p.list.3 <- list()
for(i in 1:reps){ 
  
  # randomly generate parameters, define initial state values, and set simulation length 
  ps <- c(r1 = r1 <- runif(1, min = rmin, max = rmax), 
          r2 = r2 <- runif(1, min = rmin, max = rmax), 
          K1 = K1 <- runif(1, min = Kmin, max = Kmax), 
          K2 = K2 <- runif(1, min = Kmin, max = Kmax), 
          e11 = e11 <- runif(1, min = emin, max = emax), 
          e12 = e12 <- runif(1, min = emin, max = emax), 
          e21 = e21 <- 0,#runif(1, min = emin, max = emax), 
          e22 = e22 <- 0,#runif(1, min = emin, max = emax),  
          h11 = h11 <- runif(1, min = hmin, max = hmax),
          h12 = h12 <- runif(1, min = hmin, max = hmax), 
          h21 = h21 <- 0,#runif(1, min = hmin, max = hmax), 
          h22 = h22 <- 0,#runif(1, min = hmin, max = hmax), 
          a11 = a11 <- runif(1, min = amin, max = amax), 
          a12 = a12 <- runif(1, min = amin, max = amax), 
          a21 = a21 <- 0,#runif(1, min = amin, max = amax), 
          a22 = a22 <- 0,#runif(1, min = amin, max = amax),  
          m1 = m1 <- runif(1, min = mmin, max = mmax),
          m2 = m2 <- 0,#runif(1, min = mmin, max = mmax), 
          w11 = w11 <- runif(1, min = wmin, max = wmax), 
          w22 = w22 <- 0)#runif(1, min = wmin, max = wmax))
  i.state = c(R1 = K1, R2 = K2, C1 = C.init, C2 = 0)
  
  # run simulation and create output with parameters and final state variables.
  out <- ode(parms = ps, y = i.state, times = 1:sim.length, func = ECD_model, rootfun = rootfun) %>% data.frame()
  
  feasible.3 <- out %>% assert(within_bounds(1e-3, Inf, include.lower=FALSE), R1:C1, error_fun = warning) # check for positive abundances for the 3 species throughout simulation (C2 = 0).
  
  if(class(feasible.3) == 'data.frame'){
    eq.state = tail(feasible.3, 1)
    eq.jac = jac.norm.calc(data = c(ps, R1 = eq.state$R1, R2 = eq.state$R2, C1 = eq.state$C1, C2 = eq.state$C2), jac.mat = jac) # calculate Jacobian
    
    p.list.3[[i]] <- c(ps, root.time = dim(out)[1], sim.length = sim.length, R1 = eq.state$R1, R2 = eq.state$R2, C1 = eq.state$C1, C2 = eq.state$C2, eq.max.eigen = max(Re(eigen(eq.jac[1:3,1:3])$values))) # C2 should be zero
  } else(p.list.3[[i]] <- c(ps, root.time = NA, sim.length = sim.length, R1 = NA, R2 = NA, C1 = NA, C2 = NA, eq.max.eigen = NA))
}

## Turn list into a data frame ----
df.3sp <- ldply(p.list.3) %>% tbl_df()
write.csv(x = df.3sp, file = "Google Drive/ms_ecd/df.3sp.csv")
