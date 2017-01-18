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

## maximum and minimum values for parameters ----
rmax <- 5; rmin <- 1; Kmax <- 10; Kmin <- 2; emin <- 0.5; emax <- 0.9; 
hmax <- 0.6; hmin <- 0.1; amax <- 5; amin <- 1; mmax <- 1; mmin <- 0.1; 
wmax <- 0.95; wmin <- 0.05

C.init <- 1 # initial value for consumer(s)
reps <- 50 # of parameter sets to generate
sim.length <- 1000 # duration of simulations on each parameter set

LHS.4 <- LHS(model=NULL, factors=factors, N=50, q=q, q.arg=q.arg, nboot=0) # use nboot=0 if do not need partial correlations
LHS.4.ps <- get.data(LHS.4)

LHS.4.ps <- LHS.4.ps %>% mutate(eigen = NA)

## generate for-loop for 4-species community ----
p.list.4 <- list()
for(i in 1:reps){
  
  # randomly generate parameters, define initial state values, and set simulation length 
  ps <- c(r1 = r1 <- runif(1, min = rmin, max = rmax), 
          r2 = r2 <- runif(1, min = rmin, max = rmax), 
          K1 = K1 <- runif(1, min = Kmin, max = Kmax), 
          K2 = K2 <- runif(1, min = Kmin, max = Kmax), 
          e11 = e11 <- runif(1, min = emin, max = emax), 
          e12 = e12 <- runif(1, min = emin, max = emax), 
          e21 = e21 <- runif(1, min = emin, max = emax), 
          e22 = e22 <- runif(1, min = emin, max = emax),  
          h11 = h11 <- runif(1, min = hmin, max = hmax),
          h12 = h12 <- runif(1, min = hmin, max = hmax), 
          h21 = h21 <- runif(1, min = hmin, max = hmax), 
          h22 = h22 <- runif(1, min = hmin, max = hmax), 
          a11 = a11 <- runif(1, min = amin, max = amax), 
          a12 = a12 <- runif(1, min = amin, max = amax), 
          a21 = a21 <- runif(1, min = amin, max = amax), 
          a22 = a22 <- runif(1, min = amin, max = amax),  
          m1 = m1 <- runif(1, min = mmin, max = mmax),
          m2 = m2 <- runif(1, min = mmin, max = mmax), 
          w11 = w11 <- runif(1, min = wmin, max = wmax), 
          w22 = w22 <- runif(1, min = wmin, max = wmax))
  i.state = c(R1 = K1, R2 = K2, C1 = C.init, C2 = C.init)

  # run simulation and create output with parameters and final state variables.
  out <- ode(parms = ps, y = i.state, times = 1:sim.length, func = ECD_model, rootfun = rootfun) %>% data.frame()
  
  feasible.4 <- out %>% assert(within_bounds(1e-3, Inf, include.lower=FALSE), R1:C2, error_fun = warning) # check for abundances above a threshold for all 4 species throughout simulation
  
  if(class(feasible.4) == 'data.frame'){
    eq.state = tail(feasible.4, 1)
    eq.jac = jac.norm.calc(data = c(ps, R1 = eq.state$R1, R2 = eq.state$R2, C1 = eq.state$C1, C2 = eq.state$C2), jac.mat = jac) # calculate Jacobian
    
    p.list.4[[i]] <- c(ps, root.time = dim(out)[1], sim.length = sim.length, R1 = eq.state$R1, R2 = eq.state$R2, C1 = eq.state$C1, C2 = eq.state$C2, eq.max.eigen = max(Re(eigen(eq.jac[1:4,1:4])$values)))
  } else(p.list.4[[i]] <- c(ps, root.time = NA, sim.length = sim.length, R1 = NA, R2 = NA, C1 = NA, C2 = NA, eq.max.eigen = NA))
}

## Turn list into a data frame ----
df.4sp <- ldply(p.list.4) %>% tbl_df()
as.matrix(df.4sp) # for viewing
write.csv(x = df.4sp, file = "Google Drive/ms_ecd/df.4sp.csv")

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
