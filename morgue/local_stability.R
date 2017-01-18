######################################################
## Description: 
## Code author(s): Matt Barbour
## Email: barbour@zoology.ubc.ca
######################################################

## load required libraries ----
library(deSolve)
library(plyr)
library(dplyr)

## expressions for consumer-resource dynamics WITH MUTANTS ----
    
# dynamical equations
dR1.dt <- expression(r1 * R1 * (1 - R1 / K1) - C1 * ((((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11 * R1)/(1 + ((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11 * h11 * R1 + (((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12 * h12 * R2)) - C2 * (((((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21 * R1)/(1 + (((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21 * h21 * R1 + ((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22 * h22 * R2)) - mC1 * ((((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11m * R1)/(1 + ((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11m * h11 * R1 + (((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12m * h12 * R2)) - mC2 * (((((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21m * R1)/(1 + (((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21m * h21 * R1 + ((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22m * h22 * R2)))

dR2.dt <- expression(r2 * R2 * (1 - R2 / K2) - C1 * (((((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12 * R2)/(1 + ((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11 * h11 * R1 + (((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12 * h12 * R2)) - C2 * ((((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22 * R2)/(1 + (((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21 * h21 * R1 + ((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22 * h22 * R2)) - mC1 * (((((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12m * R2)/(1 + ((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11m * h11 * R1 + (((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12m * h12 * R2)) - mC2 * ((((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22m * R2)/(1 + (((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21m * h21 * R1 + ((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22m * h22 * R2)))

dC1.dt <- expression(C1 * (e11 * ((((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11 * R1)/(1 + ((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11 * h11 * R1 + (((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12 * h12 * R2)) + e12 * (((((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12 * R2)/(1 + ((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11 * h11 * R1 + (((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12 * h12 * R2)) - m1))

dC2.dt <- expression(C2 * (e21 * (((((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21 * R1)/(1 + (((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21 * h21 * R1 + ((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22 * h22 * R2)) + e22 * ((((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22 * R2)/(1 + (((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21 * h21 * R1 + ((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22 * h22 * R2)) - m2))

dmC1.dt <- expression(mC1 * (e11 * ((((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11m * R1)/(1 + ((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11m * h11 * R1 + (((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12m * h12 * R2)) + e12 * (((((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12m * R2)/(1 + ((w11 * R1)/(w11 * R1 + (1 - w11) * R2)) * a11m * h11 * R1 + (((1 - w11) * R2)/(w11 * R1 + (1 - w11) * R2)) * a12m * h12 * R2)) - m1))

dmC2.dt <- expression(mC2 * (e21 * (((((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21m * R1)/(1 + (((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21m * h21 * R1 + ((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22m * h22 * R2)) + e22 * ((((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22m * R2)/(1 + (((1 - w22) * R1)/((1 - w22) * R1 + w22 * R2)) * a21m * h21 * R1 + ((w22 * R2)/((1 - w22) * R1 + w22 * R2)) * a22m * h22 * R2)) - m2))


## setup Jacobian ----
jac <- list(D(dR1.dt,"R1"), D(dR2.dt,"R1"), D(dC1.dt,"R1"), D(dC2.dt,"R1"), D(dmC1.dt,"R1"), D(dmC2.dt,"R1"),
            D(dR1.dt,"R2"), D(dR2.dt,"R2"), D(dC1.dt,"R2"), D(dC2.dt,"R2"), D(dmC1.dt,"R2"), D(dmC2.dt,"R2"),
            D(dR1.dt,"C1"), D(dR2.dt,"C1"), D(dC1.dt,"C1"), D(dC2.dt,"C1"), D(dmC1.dt,"C1"), D(dmC2.dt,"C1"),
            D(dR1.dt,"C2"), D(dR2.dt,"C2"), D(dC1.dt,"C2"), D(dC2.dt,"C2"), D(dmC1.dt,"C2"), D(dmC2.dt,"C2"),
            D(dR1.dt,"mC1"), D(dR2.dt,"mC1"), D(dC1.dt,"mC1"), D(dC2.dt,"mC1"), D(dmC1.dt,"mC1"), D(dmC2.dt,"mC1"),
            D(dR1.dt,"mC2"), D(dR2.dt,"mC2"), D(dC1.dt,"mC2"), D(dC2.dt,"mC2"), D(dmC1.dt,"mC2"), D(dmC2.dt,"mC2"))

## compute jacobian for invasion of mutant ----
jac.norm.calc <- function(data, jac.mat){
  
  # assign parameters and equilibrium state variables to numerically calculate the Jacobian
  r1 <- data["r1"]
  r2 <- data["r2"]
  K1 <- data["K1"]
  K2 <- data["K2"]
  e11 <- data["e11"]
  e12 <- data["e12"]
  e21 <- data["e21"]
  e22 <- data["e22"]
  h11 <- data["h11"]
  h12 <- data["h12"]
  h21 <- data["h21"]
  h22 <- data["h22"]
  a11 <- data["a11"]
  a12 <- data["a12"]
  a21 <- data["a21"]
  a22 <- data["a22"]
  a11m <- 0
  a12m <- 0
  a21m <- 0
  a22m <- 0
  m1 <- data["m1"]
  m2 <- data["m2"]
  w11 <- data["w11"]
  w22 <- data["w22"]
  R1 <- data["R1"]
  R2 <- data["R2"]
  C1 <- data["C1"]
  C2 <- data["C2"]
  mC1 <- 0 
  mC2 <- 0
  
  jac.norm.calc <- matrix(sapply(jac.mat, function(pd) eval(pd)), nrow = 6)
}

## calculate maximum Real eigenvalue
max.eigen <- function(data, jac.mat){
  eq.jac <- jac.norm.calc(data, jac.mat)
  eq.jac.mat <- matrix(unlist(eq.jac), nrow = 6)
  return(max(Re(eigen(eq.jac.mat[1:4,1:4])$values)))
}

## compute jacobian for invasion of mutant ----
jac.mut.calc <- function(data, jac.mat){
  
  # assign parameters and equilibrium state variables to numerically calculate the Jacobian
  r1 <- data["r1"]
  r2 <- data["r2"]
  K1 <- data["K1"]
  K2 <- data["K2"]
  e11 <- data["e11"]
  e12 <- data["e12"]
  e21 <- data["e21"]
  e22 <- data["e22"]
  h11 <- data["h11"]
  h12 <- data["h12"]
  h21 <- data["h21"]
  h22 <- data["h22"]
  a11 <- data["a11"]
  a12 <- data["a12"]
  a21 <- data["a21"]
  a22 <- data["a22"]
  a11m <- data["a11m"]
  a12m <- data["a12m"]
  a21m <- data["a21m"]
  a22m <- data["a22m"]
  m1 <- data["m1"]
  m2 <- data["m2"]
  w11 <- data["w11"]
  w22 <- data["w22"]
  R1 <- data["R1"]
  R2 <- data["R2"]
  C1 <- data["C1"]
  C2 <- data["C2"]
  mC1 <- 0 
  mC2 <- 0
  
  jac.mut.calc <- matrix(sapply(jac.mat, function(pd) eval(pd)), nrow = 6)
}
