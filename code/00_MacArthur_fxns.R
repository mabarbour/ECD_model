######################################################
## Description: Dynamical model and other useful functions for checking dynamics
## Code author(s): Matt Barbour
## Email: barbour@zoology.ubc.ca
######################################################


## Consumer-resource dynamics ----
# general assumptions:
# (1) logistic growth for the resource(s)
# (2) multispecies, spatially-implicit type 2 functional response of consumer 
# (3) density-independent mortality for consumer 
# function is in the necessary format "function(Time,State,Pars)" for solving used the ode() function in the R package "deSolve"

## 4 species model: 2 consumers, 2 resources
ECD_model_Mac <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    # consumer functional responses
    C1R1fxn <- a11 * R1 # (W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- a12 * R2 # (W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C2R1fxn <- a21 * R1 # (W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    C2R2fxn <- a22 * R2 # (W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - C2 * C2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - C2 * C2R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## 3 species model: C1, R1, R2
ECD_model.C1.3sp_Mac <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    
    # consumer functional responses
    C1R1fxn <- a11 * R1 #(W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- a12 * R2 #(W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    #C2R1fxn <- (W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    #C2R2fxn <- (W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn #- C2 * C2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn #- C2 * C2R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    #dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt))) #, dC2.dt
  })
}

## 3 species model: C2, R1, R2
ECD_model.C2.3sp_Mac <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    
    # consumer functional responses
    #C1R1fxn <- (W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    #C1R2fxn <- (W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C2R1fxn <- a21 * R1 #(W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    C2R2fxn <- a22 * R2 #(W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C2 * C2R1fxn #- C1 * C1R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C2 * C2R2fxn #- C1 * C1R2fxn
    #dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC2.dt))) #, dC1.dt
  })
}


## With mutants - Consumer-resource dynamics ----
# same assumptions as for ECD_model

## 4 species model with possible C1 and C2 mutants
mut_ECD_model_Mac <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    # normal consumer functional responses
    C1R1fxn <- a11 * R1 #(W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- a12 * R2 #(W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C2R1fxn <- a21 * R1 #(W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    C2R2fxn <- a22 * R2 #(W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # mutant consumer functional responses
    mC1R1fxn <- a11m * R1 #(W11 * a11m * R1)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    mC1R2fxn <- a12m * R2 #(W12 * a12m * R2)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    mC2R1fxn <- a21m * R1 #(W21 * a21m * R1)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    mC2R2fxn <- a22m * R2 #(W22 * a22m * R2)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - C2 * C2R1fxn - mC1 * mC1R1fxn - mC2 * mC2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - C2 * C2R2fxn - mC1 * mC1R2fxn - mC2 * mC2R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    dmC1.dt <- mC1 * (e11 * mC1R1fxn + e12 * mC1R2fxn - m1)
    dmC2.dt <- mC2 * (e21 * mC2R1fxn + e22 * mC2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt, dmC1.dt, dmC2.dt)))
  })
}

## 4 species model with only C1 mutant
mutC1_ECD_model_Mac <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    # normal consumer functional responses
    C1R1fxn <- a11 * R1 #(W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- a12 * R2 #(W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C2R1fxn <- a21 * R1 #(W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    C2R2fxn <- a22 * R2 #(W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # mutant consumer functional responses
    mC1R1fxn <- a11m * R1 #(W11 * a11m * R1)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    mC1R2fxn <- a12m * R2 #(W12 * a12m * R2)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    #mC2R1fxn <- (W21 * a21m * R1)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    #mC2R2fxn <- (W22 * a22m * R2)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - C2 * C2R1fxn - mC1 * mC1R1fxn #- mC2 * mC2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - C2 * C2R2fxn - mC1 * mC1R2fxn #- mC2 * mC2R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    dmC1.dt <- mC1 * (e11 * mC1R1fxn + e12 * mC1R2fxn - m1)
    #dmC2.dt <- mC2 * (e21 * mC2R1fxn + e22 * mC2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt, dmC1.dt))) # dmC2.dt
  })
}

## 4 species model with only C2 mutant
mutC2_ECD_model_Mac <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    # normal consumer functional responses
    C1R1fxn <- a11 * R1 #(W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- a12 * R2 #(W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C2R1fxn <- a21 * R1 #(W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    C2R2fxn <- a22 * R2 #(W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # mutant consumer functional responses
    #mC1R1fxn <- (W11 * a11m * R1)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    #mC1R2fxn <- (W12 * a12m * R2)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    mC2R1fxn <- a21m * R1 #(W21 * a21m * R1)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    mC2R2fxn <- a22m * R2 #(W22 * a22m * R2)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - C2 * C2R1fxn - mC2 * mC2R1fxn # - mC1 * mC1R1fxn 
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - C2 * C2R2fxn - mC2 * mC2R2fxn # - mC1 * mC1R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    #dmC1.dt <- mC1 * (e11 * mC1R1fxn + e12 * mC1R2fxn - m1)
    dmC2.dt <- mC2 * (e21 * mC2R1fxn + e22 * mC2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt, dmC2.dt))) # dmC1.dt
  })
}

## 3 species model with only C1 mutant
mut_ECD_model.C1.3sp_Mac <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    
    # normal consumer functional responses
    C1R1fxn <- a11 * R1 #(W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- a12 * R2 #(W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    #C2R1fxn <- (W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    #C2R2fxn <- (W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # mutant consumer functional responses
    mC1R1fxn <- a11m * R1 #(W11 * a11m * R1)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    mC1R2fxn <- a12m * R2 #(W12 * a12m * R2)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    #mC2R1fxn <- (W21 * a21m * R1)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    #mC2R2fxn <- (W22 * a22m * R2)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - mC1 * mC1R1fxn # - C2 * C2R1fxn - mC2 * mC2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - mC1 * mC1R2fxn # - C2 * C2R2fxn - mC2 * mC2R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    #dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    dmC1.dt <- mC1 * (e11 * mC1R1fxn + e12 * mC1R2fxn - m1)
    #dmC2.dt <- mC2 * (e21 * mC2R1fxn + e22 * mC2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dmC1.dt))) # dC2.dt, dmC2.dt
  })
}

## 3 species model with only C2 mutant
mut_ECD_model.C2.3sp_Mac <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    
    # normal consumer functional responses
    #C1R1fxn <- (W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    #C1R2fxn <- (W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C2R1fxn <- a21 * R1 #(W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    C2R2fxn <- a22 * R2 #(W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # mutant consumer functional responses
    #mC1R1fxn <- (W11 * a11m * R1)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    #mC1R2fxn <- (W12 * a12m * R2)/(1 + W11 * a11m * h11 * R1 + W12 * a12m * h12 * R2)
    mC2R1fxn <- a21m * R1 #(W21 * a21m * R1)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    mC2R2fxn <- a22m * R2 #(W22 * a22m * R2)/(1 + W21 * a21m * h21 * R1 + W22 * a22m * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C2 * C2R1fxn - mC2 * mC2R1fxn # - C1 * C1R1fxn - mC1 * mC1R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C2 * C2R2fxn - mC2 * mC2R2fxn # - C1 * C1R2fxn - mC1 * mC1R2fxn 
    #dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    #dmC1.dt <- mC1 * (e11 * mC1R1fxn + e12 * mC1R2fxn - m1)
    dmC2.dt <- mC2 * (e21 * mC2R1fxn + e22 * mC2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC2.dt, dmC2.dt))) # dC1.dt, dmC1.dt
  })
}


## Identify steady-state equilibriums ----
# using failwith() to make function give output even with a fatal error
# runsteady() identifies the steady-state of the model.
safe.runsteady <- failwith(default = structure(list(y = c(R1 = NA, R2 = NA, C1 = NA, C2 = NA)), steady = 0), runsteady.ps)
#safe.runsteady.C1 <- failwith(default = c(C1.R1 = NA, C1.R2 = NA, C1.C1 = NA, steady = 0), runsteady.ps)
#safe.runsteady.C2 <- failwith(default = c(C2.R1 = NA, C2.R2 = NA, C2.C2 = NA, steady = 0), runsteady.ps)
#safe.jacobian <- failwith()


## Plotting consumer-resource dynamics over time ----
dynamic_matplot <- function(param.vector, # vector of parameters for the dynamical model
                            init.state, # initial state variables for the model
                            sim.length, # length of simulation
                            model, # dynamical model
                            ylim, # limits of state variables for plotting
                            ...){
  
  # Run the experiment
  df <- ode(init.state, 1:sim.length, model, param.vector)
  
  # plot the results. 
  matplot(df[ ,"time"], df[ ,names(init.state)],
          type = "l", ylim=ylim, ...)
}
