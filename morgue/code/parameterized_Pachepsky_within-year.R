
##############################################################################################################
## Estimate stickleback consumer-resource dynamics based on allometric scaling and literature.
## Author: Matt Barbour
## Email: matthew.a.barbour@gmail.com
##############################################################################################################


## LOAD REQUIRE FUNCTIONS ----

source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')
source('code/00_Pachepsky_fxns.R')


mass_scaling <- function(W, intercept, exponent) {
  log_rate <- intercept + exponent * log(W)
  return(exp(log_rate))
}


## BODY SIZE (g) ----

g_zoop <- 0.004/1000            # covert 0.004 mg to g
g_benthic_invert <- 0.25/1000   # convert 0.25 mg to g

g_limnetic <- 0.81              # from Schluter 1993, Ecology Table 1
g_benthic <- 1.48               # from Schluter 1993, Ecology Table 1
g_solitary <- 1.12              # from Schluter 1993, Ecology Table 1


## INTRINSIC GROWTH RATE (r) ----
# r = birth rate per day
# proportional to M^(-1/4) from Rip and McCann 2011 (they cite Brown et al. 2004)
# alternatively use scaling exponent = -0.27377 for heterotherm metazoa from Fenchel 1974, Oecologia
r_exponent <- -0.25             # theoretical scaling.
r_scale <- 0.3                  # estimate of r from Daphnia pulex in Table 1 of Fenchel 1974. Original study is Frank et al. 1957       
r_intercept <- log(r_scale) - r_exponent * log(g_zoop) # scale based on zooplankton size

r_benthos_estimate <- mass_scaling(W = g_benthic_invert, intercept = r_intercept, exponent = r_exponent) 


#### CARRYING CAPACITY (K) ----
# proportional to M^(-3/4) from Rip and McCann 2011 (they cite Enquist et al. 1998 and Brown et al. 2004)
# alternatively use empirical exponent = -0.77 from Brown et al. 2004
K_exponent <- -0.75 
K_scale <- 5          # admittedly, this is a bit of an arbitrary number. It does correspond to zooplankton biomass density (g/500 L) from fishless ponds
K_intercept <- log(K_scale) - K_exponent * log(g_zoop) # scale based on zooplankton size

K_benthos_estimate <- mass_scaling(W = g_benthic_invert, intercept = K_intercept, exponent = K_exponent)


#### CONVERSION EFFICIENCY (e) ----
# from Mittelbach 1981 (see refs to Ware 1975 and Elliot 1976) 
# independent of body size (Rip and McCann 2011, Yodzis & Innes 1992)
e <- 0.7 


#### HABITAT PREFERENCE (w) ----
# Based on % littoral Carbon from Matthews et al. 2010
limnetics_littoral <- c(0.4, 0.1)
benthics_littoral <- c(0.7, 0.6)
solitary_littoral <- c(0.5, 0.6, 1.0, 0.1)

w_limnetics_in_open_water <- 1 - mean(limnetics_littoral) # 0.65 is approximated from Fig. 2 in Schluter 1993
w_benthics_in_benthos <- mean(benthics_littoral) # 0.9 is approximated from Fig. 2 in Schluter 1993       

# as 1 - littoral to correspond to parameters in ECD_model.C1.3sp 
w_solitary <- 1 - mean(solitary_littoral)                 # 0.5

#### MORTALITY RATE (m) ----
# m = individual deaths per unit time
# proportional to M^(-1/4) from Brown et al. 2004
# empirical estimate of exponent = -0.24 from Fig. 4 in Brown et al. 2004

# estimate of annual threespine stickleback mortality from fishbase
library(rfishbase) 
threespine_stickleback_M <- filter(popgrowth(species_list = "Gasterosteus aculeatus"), M != "NA")$M 

m_exponent <- -0.25                   # theoretical exponent
m_scale <- threespine_stickleback_M/365        
m_intercept <- log(m_scale) - m_exponent * log(g_solitary) # use solitary stickleback mass for scaling

m_limnetic_estimate <- mass_scaling(W = g_limnetic, intercept = m_intercept, exponent = m_exponent); m_limnetic_estimate
m_benthic_estimate <- mass_scaling(W = g_benthic, intercept =  m_intercept, exponent = m_exponent); m_benthic_estimate


#### HANDLING TIME (h) ---

# from Mittelbach 1981, Table 3: Minimum handling time for Daphnia vs. Chironomus. Units = seconds per individual
min_h_Chironomus <- 9.63  # this can be much higher because stickleback than sunfish
min_h_Daphnia <- 1.02
benthos_relative_handling_time <- min_h_Chironomus/min_h_Daphnia 

h_scale <- 0.01 # assuming it's independent of consumer body size. This assumption is likely okay, given that prey varied 2-orders of magnitude while stickleback vary much less than 1-order.


#### ATTACK RATE (a) ----
# estimates from Fig. 4, Schluter 1993
# although these are estimates of prey volume consumed per strike, 
# Schluter noted that these relative estimates were virtually identical to intake rate.
# The intake rate data taken by Schluter  is more in line with per capita attack rate; however,
# not all of these data were reported. This is why we are using capture success as a proxy.
# We use these estimates instead of mass-based estimates of attack rate (e.g. Rip and McCann 2011 cite Peters 1983 and Brown et al. 2004),
# since these traits depend other evolved traits that are not related to body size (e.g. gill rakers, gape size and body shape)

limnetic_capture_success_open_water <- 0.011
limnetic_capture_success_benthos <- 0.025        # note that limnetics have an overall higher capture succes in benthos vs. open water, this doesn't nullify the existence of a trade-off between the two habitats
benthic_capture_success_open_water <- 0.004 
benthic_capture_success_benthos <- 0.16

solitary_capture_success_open_water <- 0.009
solitary_capture_success_benthos <- 0.04

a_scale <- 50#300                                   # scaled until I was able to get coexistence for species pairs 


#### CONSUMER-RESOURCE DYNAMICS ----

## Species-pair lakes
species_pair.init.states <- c(R1=1, R2=1, C1=0.25, C2=0.25) # note that equilibrium dynamics appear to be independent of initial state variables
species_pair.parms <- data.frame(r1 = r_scale, 
                                 r2 = r_benthos_estimate, 
                                 K1 = K_scale, 
                                 K2 = K_benthos_estimate,
                                 e11 = e, 
                                 e12 = e, 
                                 e21 = e, 
                                 e22 = e,
                                 h11 = h_scale, 
                                 h12 = h_scale*benthos_relative_handling_time, 
                                 h21 = h_scale, 
                                 h22 = h_scale*benthos_relative_handling_time,
                                 a11 = a_scale*limnetic_capture_success_open_water, 
                                 a12 = a_scale*limnetic_capture_success_benthos, 
                                 a21 = a_scale*benthic_capture_success_open_water, 
                                 a22 = a_scale*benthic_capture_success_benthos,
                                 w11 = w_limnetics_in_open_water, 
                                 w22 = w_benthics_in_benthos, 
                                 m1 = m_limnetic_estimate, 
                                 m2 = m_benthic_estimate)

species_pair_dynamics <- safe.runsteady(y = species_pair.init.states, func = Pachepsky_model, parms = species_pair.parms, stol = 1e-5)
species_pair_dynamics$y # equilibrium densities

# plot the dynamics
matplot.0D(ode(species_pair.init.states, 1:500, Pachepsky_model, species_pair.parms), ylim=c(0,K_scale))


## Solitary lakes
solitary.init.states <- c(R1=1, R2=1, C1=0.5) # note that equilibrium dynamics appear to be independent of initial state variables
solitary.parms <- data.frame(r1 = r_scale, 
                             r2 = r_benthos_estimate, 
                             K1 = K_scale, 
                             K2 = K_benthos_estimate,
                             e11 = e, 
                             e12 = e, 
                             #e21 = e, 
                             #e22 = e,
                             h11 = h_scale, 
                             h12 = h_scale*benthos_relative_handling_time, 
                             #h21 = h_scale, 
                             #h22 = h_scale*benthos_relative_handling_time,
                             a11 = a_scale*solitary_capture_success_open_water, 
                             a12 = a_scale*solitary_capture_success_benthos, 
                             #a21 = a_scale*benthic_capture_success_open_water, 
                             #a22 = a_scale*benthic_capture_success_benthos,
                             w11 = w_solitary, 
                             #w22 = w_benthics_in_benthos, 
                             m1 = m_scale) 
#m2 = m_benthic_estimate)

solitary_dynamics <- safe.runsteady(y = solitary.init.states, func = Pachepsky_model.C1.3sp, parms = solitary.parms, stol = 1e-5)
solitary_dynamics$y # equilibrium densities

# plot the dynamics
matplot.0D(ode(solitary.init.states, 1:500, Pachepsky_model.C1.3sp, solitary.parms), ylim=c(0,K_scale))


## COMPARE SPECIES-PAIR & SOLITARY LAKES ---- 

## Resources

# relative densities of total resources
(species_pair_dynamics$y["R1"] + species_pair_dynamics$y["R2"])/(solitary_dynamics$y["R1"] + solitary_dynamics$y["R2"]) - 1 

# relative densities of zooplankton
species_pair_dynamics$y["R1"]/solitary_dynamics$y["R1"] - 1 # zooplankton densities are lower in species pair lakes

# relative densities of benthic inverts
species_pair_dynamics$y["R2"]/solitary_dynamics$y["R2"] - 1 

## Consumers

# relative densities of total consumers
(species_pair_dynamics$y["C1"] + species_pair_dynamics$y["C2"])/solitary_dynamics$y["C1"] - 1

# relative densities of limnetics compare to benthics
species_pair_dynamics$y["C1"]/species_pair_dynamics$y["C2"]

## Consumer:Resource Ratios
log(solitary_dynamics$y["C1"]/((solitary_dynamics$y["R1"] + solitary_dynamics$y["R2"])))
log((species_pair_dynamics$y["C1"] + species_pair_dynamics$y["C2"])/((species_pair_dynamics$y["R1"] + species_pair_dynamics$y["R2"])))


