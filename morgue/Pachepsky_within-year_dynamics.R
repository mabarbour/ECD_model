

## source general parameters and functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')
source('code/00_Pachepsky_fxns.R')

## general parameters for simulations
r <- 1
K <- 4
#e <- 0.8
h <- 0.4
aii <- 10
aij <- 1
w <- 0.8
m <- 0.001
R <- 2
C <- 1

## Species-pair lakes
species_pair.init.states <- c(R1=1, R2=1, C1=0.25, C2=0.25) # note that equilibrium dynamics appear to be independent of initial state variables
species_pair.parms <- data.frame(r1 = r, r2 = r, K1 = K, K2 = K,
                                 #e11 = e, e12 = e, e21 = e, e22 = e,
                                 h11 = h, h12 = h, h21 = h, h22 = h,
                                 a11 = aii, a12 = aij, a21 = aij, a22 = aii,
                                 w11 = w, w22 = w, m1 = m, m2 = m)

species_pair_dynamics <- safe.runsteady(y = species_pair.init.states, func = Pachepsky_model, parms = species_pair.parms, stol = 1e-5)
species_pair_dynamics$y # equilibrium densities

# plot the dynamics
matplot.0D(ode(species_pair.init.states, 1:100, Pachepsky_model, species_pair.parms), ylim=c(0,K))


## Solitary lakes
w_solitary <- 0.5
aii_solitary <- (aii + aij)/2
aij_solitary <- aii_solitary

solitary.init.states <- c(R1=1, R2=1, C1=0.5) # note that equilibrium dynamics appear to be independent of initial state variables
solitary.parms <- data.frame(r1 = r, r2 = r, K1 = K, K2 = K,
                             #e11 = e, e12 = e, e21 = e, e22 = e,
                             h11 = h, h12 = h,
                             a11 = aii_solitary, a12 = aij_solitary,
                             w11 = w_solitary, m1 = m) 

solitary_dynamics <- safe.runsteady(y = solitary.init.states, func = Pachepsky_model.C1.3sp, parms = solitary.parms, stol = 1e-5)
solitary_dynamics$y # equilibrium densities

# plot the dynamics
matplot.0D(ode(solitary.init.states, 1:100, Pachepsky_model.C1.3sp, solitary.parms), ylim=c(0,K)) #K_scale))


## COMPARE SPECIES-PAIR & SOLITARY LAKES ---- 

## Resources

# relative densities of total resources
(species_pair_dynamics$y["R1"] + species_pair_dynamics$y["R2"])/(solitary_dynamics$y["R1"] + solitary_dynamics$y["R2"]) - 1 

# relative densites of zooplankton
species_pair_dynamics$y["R1"]/solitary_dynamics$y["R1"] - 1 

# relative densities of benthic inverts
species_pair_dynamics$y["R2"]/solitary_dynamics$y["R2"] - 1 

## Consumers

# relative densities of total consumers 
(species_pair_dynamics$y["C1"] + species_pair_dynamics$y["C2"])/solitary_dynamics$y["C1"] - 1

# relative densities of limentics to benthics
species_pair_dynamics$y["C1"]/species_pair_dynamics$y["C2"]

## Consumer:Resource Ratios
log(solitary_dynamics$y["C1"]/((solitary_dynamics$y["R1"] + solitary_dynamics$y["R2"])))
log((species_pair_dynamics$y["C1"] + species_pair_dynamics$y["C2"])/((species_pair_dynamics$y["R1"] + species_pair_dynamics$y["R2"])))



