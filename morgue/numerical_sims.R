######################################################
## Description: 
## Code author(s): Matt Barbour
## Email: barbour@zoology.ubc.ca
######################################################

## load required libraries ----
library(deSolve)
library(primer)
library(plyr)
library(dplyr)
library(assertr)

## source required functions ----
source('Google Drive/ms_ecd/numerical_fxns.R')
source('Google Drive/ms_ecd/local_stability.R')

## load data
df.4sp <- read.csv("Google Drive/ms_ecd/hypercube.4sp.csv") %>% tbl_df()
df.4sp.feas <- df.4sp %>% select(-X) %>% filter(R1>1e-3, R2>1e-3, C1>1e-3, C2>1e-3, steady == 1) # keep feasible and stable equilibriums

## choose parameter set to examine
set <- as.matrix(df.4sp.feas[1, ])[1, ]

## plot dynamics
#dynamic_matplot(param.vector = set[5:24],
 #               init.state = c(R1 = R1 <- set[1], 
  #                             R2 = R2 <- set[2], 
   #                            C1 = C1 <- set[3], 
    #                           C2 = C2 <- set[4]),
     #           sim.length = 1000,
      #          model = ECD_model, 
       #         ylim = c(0,max(set[c("K1","K2")])))

## setup simulation
sim.length <- 1000
reps <- 50 # number of mutation iterations
sim.list.4 <- list() 
sim.list.4[[1]] <- c(set, eq.max.eigen = NA) # input first data vector

## run simulation
for(i in 2:reps){
  
  # extract parameters and state variables
  ps <- sim.list.4[[i-1]][1:20]
  i.state <- sim.list.4[[i-1]][23:26]
  
  # input parameter values and state variable equilibriums and run the simulation
  out1 <- ode(parms = ps, y = i.state, times = 1:sim.length, func = ECD_model, rootfun = rootfun) %>% data.frame()
  
  # assign mutation in attack rate to one of the species
  sp <- rbinom(1,1,0.5)
  mut <- ifelse(rbinom(1,1,0.5) == 0, -0.1, 0.1)
  
  # determine whether mutant can invade or not by calculating whether it has positive growth rate (do based on simplified expression...or based on eigenvalue using the matrix approach?)
  new.ps <- ps
  
  if(sp == 0){
    new.ps["a11m"] <- new.ps["a11"] + mut
    new.ps["a12m"] <- new.ps["a12"] - mut
    new.ps["a21m"] <- 0
    new.ps["a22m"] <- 0
    jac.mut <- jac.mut.calc(data = c(new.ps, i.state), jac.mat = jac)
    eigen.mut <- Re(eigen(jac.mut[5,5])$values)
    #max.int.eigen <- max(Re(eigen(jac.mut[1:4,1:4])$values))
  }
  if(sp == 1){
    new.ps["a11m"] <- 0
    new.ps["a12m"] <- 0
    new.ps["a22m"] <- new.ps["a22"] + mut
    new.ps["a21m"] <- new.ps["a21"] - mut
    jac.mut <- jac.mut.calc(data = c(new.ps, i.state), jac.mat = jac)
    eigen.mut <- Re(eigen(jac.mut[6,6])$values)
    #max.int.eigen <- max(Re(eigen(jac.mut[1:4,1:4])$values))
  }

  # if the mutant can invade, update attack rates and print new equilibrium values
  if(eigen.mut > 0){
    if(sp == 0){
      ps["a11"] <- new.ps["a11m"]
      ps["a12"] <- new.ps["a12m"]
    }
    if(sp == 1){
      ps["a21"] <- new.ps["a21m"]
      ps["a22"] <- new.ps["a22m"]
    }
    out2 <- ode(parms = ps, y = i.state, times = 1:sim.length, func = ECD_model, rootfun = rootfun) %>% 
      data.frame() 
  }
  # if the mutant can't invade, keep attack rates the same and print the same equilibrium values.
  if(eigen.mut < 0){
    out2 <- ode(parms = ps, y = i.state, times = 1:sim.length, func = ECD_model, rootfun = rootfun) %>% 
      data.frame()
  }

  # check that new equilibrium is still feasible
  feas <- out2 %>% assert(within_bounds(0, Inf, include.lower=FALSE), R1:C2, error_fun = warning) 
  
  if(class(feas) == 'data.frame'){
    eq.state = tail(feas, 1)
    eq.jac = jac.norm.calc(data = c(ps, R1 = eq.state$R1, R2 = eq.state$R2, C1 = eq.state$C1, C2 = eq.state$C2), jac.mat = jac) # calculate Jacobian
    
    sim.list.4[[i]] <- c(ps, root.time = dim(feas)[1], sim.length = sim.length, R1 = eq.state$R1, R2 = eq.state$R2, C1 = eq.state$C1, C2 = eq.state$C2, eq.max.eigen = max(Re(eigen(eq.jac[1:4,1:4])$values))) #which.eq.max.eigen = which.max(Re(eigen(eq.jac[1:4,1:4])$values))
  } else(sim.list.4[[i]] <- c(ps, root.time = NA, sim.length = sim.length, R1 = NA, R2 = NA, C1 = NA, C2 = NA, eq.max.eigen = NA))
}

sim.4.df <- sim.list.4 %>% ldply() %>% mutate(reps = 1:reps) 

sim.4.df %>% assert(within_bounds(0, Inf, include.lower=FALSE), r1:C2, error_fun = warning) # check for positive parameter values and species abundances. 

## plots over feasible range

# stability
plot(-1*eq.max.eigen ~ reps, sim.4.df[1:39,]) # stability decreases over long term (goes closer to zero)
scatter.smooth(-1*sim.4.df[1:39,]$eq.max.eigen, sim.4.df[1:39,]$reps)

# attack rate
# interesting that both consumers evolve higher attack rates on resource 2, this is because it has a much higher reproductive rate.
plot(a12/(a11 + a12) ~ reps, sim.4.df[1:32, ]) # increase in specialization on resource 2
plot(a22/(a21 + a22) ~ reps, sim.4.df[1:32, ]) # increase in specialization on resource 2


## What do I need?
# Give a parameter set to the model, calculate the equilibrium values for resources and consumers. Then, calculate the local stability of that model.
# also, given the equilibrium values and internal stability, I would like to see if a mutant consumer (C1 or C2) could invade the system. Therefore, for a given initial parameter set, we can see how evolution would proceed for that system. We can then evaluate whether evolution, via trait divergence, eventually leads to instability in the system.