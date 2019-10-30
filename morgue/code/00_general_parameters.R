
## load required libraries
library(deSolve) # numerical integration library
library(rootSolve) # for runsteady integration function
library(pse) # for latin hypercube
library(plyr) # for ldply
library(dplyr) # for manipulating data
library(tidyr) # for tidying data
library(cowplot) # for better base, ggplot graphics

## Abundance threshold for feasible equilibriums
#abund.thres <- 1e-2 

## Trait change threshold for identifying ESS
#trait.thres <- 1e-4

## Number of mutation iterations
#mut.reps <- 10000 

## Mutation size
#mut.size <- 0.01

## All runsteady simulations go for 1000 time steps with a steady state tolerance of 1e-4. 
runsteady.ps <- function(y, func = func, parms = parms, times = c(0,1000), stol = 1e-4){
  runsteady(y = y, times = times, func = func, parms = parms, stol = stol)
}