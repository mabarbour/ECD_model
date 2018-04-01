
## load required functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

library(rfishbase)

mass_scaling <- function(W, intercept, exponent) {
  log_rate <- intercept + exponent * log(W)
  return(exp(log_rate))
}

## BODY SIZE (g) ----

g_Daphnia_pulex <- 1.6 * 10e-4  # from Table 1 in Fenchel 1974 (they cite Frank et al. 1957)
g_Chironomid <- 1.6 * 10e-3 # very rough estimate, need to confirm

g_limnetic <- 0.81 # from Schluter 1993, Ecology Table 1
g_benthic <- 1.48  # from Schluter 1993, Ecology Table 1
g_solitary <- 1.12 # from Schluter 1993, Ecology Table 1


## INTRINSIC GROWTH RATE (r) ----
# r = birth rate per day
# proportional to M^(-1/4) from Rip and McCann 2011 (they cite Brown et al. 2004)
# using relationship for heterotherm metazoa from Fenchel 1974, Oecologia

r_Daphnia_pulex <- 0.3          # from Table 1 in Fenchel 1974 (they cite Frank et al. 1957)

r_zoops_estimate <- mass_scaling(W = g_Daphnia_pulex, intercept = -1.6391, exponent = -0.2738) 
r_zoops_estimate/r_Daphnia_pulex # much higher than Frank et al. 1957 estimate...I'm going to stick with the mass-based estimates though to get the relative differences accurate between zooplankton and benthic inverts

r_benthos_estimate <- mass_scaling(W = g_Chironomid, intercept = -1.6391, exponent = -0.2738) 


#### CARRYING CAPACITY (K) ----
# proportional to M^(-3/4) from Rip and McCann 2011 (they cite Enquist et al. 1998 and Brown et al. 2004)

K_intercept <- -14.29 # from Brown et al. 2004, Fig. 6
K_zoops_estimate <- mass_scaling(W = g_Daphnia_pulex, intercept = K_intercept, exponent = -0.77); K_zoops_estimate
K_benthos_estimate <- mass_scaling(W = g_Chironomid, intercept = K_intercept, exponent = -0.77); K_benthos_estimate


#### CONVERSION EFFICIENCY (e) ----
# from Mittelbach 1981 (see refs to Ware 1975 and Elliot 1976) 
# independent of body size (Rip and McCann 2011, Yodzis & Innes 1992)
e <- 0.7 


#### HABITAT PREFERENCE (w) ----
# Based on % littoral Carbon from Matthews et al.
w_limnetics_in_open_water <- 0.80
w_benthics_in_benthos <- 0.75

w_solitary <- 0.5

#### MORTALITY RATE (m) ----
# temperature-corrected mortality rate, ln(Ze^E/kT) measured per year  and body mass, ln(M) measured in grams
E <- 0.63 # units = eV
k <- 8.617 * 10e-5 # eV/K
Celsius_to_Kelvin <- function(C) C + 273.15
-E/k*Celsius_to_Kelvin(19)
exp(1)^(E/k*Celsius_to_Kelvin(19))

temp_adj_mortality <- function(W){
  y <- -0.24*log(W) + 26.25
  return(exp(y))
}
temp_adj_mortality(g_limnetic)/365

m_intercept <- 0.4 # arbitrarily chosen
m_limnetic_estimate <- mass_scaling(W = g_limnetic, intercept = m_intercept, exponent = -0.25); m_limnetic_estimate
m_benthic_estimate <- mass_scaling(W = g_benthic, intercept =  m_intercept, exponent = -0.25); m_benthic_estimate

m_solitary_estimate <- mass_scaling(W = g_solitary, intercept = m_intercept, exponent = -0.25); m_solitary_estimate


#### HANDLING TIME (h) ---
# from Mittelbach 1981, Table 3: Minimum handling time for Daphnia vs. Chironomus
min_h_Chironomus <- 9.63 # seconds
min_h_Daphnia <- 1.02    # seconds
benthos_relative_handling_time <- min_h_Chironomus/min_h_Daphnia 

h_scale <- 0.01


#### ATTACK RATE (a) ----
# estimates from Fig. 4, Schluter 1993
# although these are estimates of prey volume consumed per strike, 
# Schluter noted that these relative estimates were virtually identical to intake rate.
# The intake rate data taken by Schluter  is more in line with per capita attack rate; however,
# not all of these data were reported. This is why we are using capture success as a proxy.
# We use these estimates instead of mass-based estimates of attack rate (e.g. Rip and McCann 2011 cite Peters 1983 and Brown et al. 2004),
# since these traits depend other evolved traits that are not related to body size (e.g. gill rakers, ...)

# note that limnetics have an overall higher capture succes in benthos vs. open water, this doesn't nullify the 
# existence of a trade-off between the two habitats

# alternatively...use Seth's estimates in g/day
limnetic_capture_success_open_water <- 0.35 # 0.011 
limnetic_capture_success_benthos <- 0.20 # 0.025
benthic_capture_success_open_water <- 0.129
benthic_capture_success_benthos <- 0.73

solitary_capture_success_open_water <- 0.181
solitary_capture_success_benthos <- 0.34

a_scale <- 2 # 100

#### DEFINE INITIAL STATE VARIABLES & PARAMETERS ----

## Species-pair food web
species_pair.init.states <- c(R1=1, R2=1, C1=0.25, C2=0.25) # note that equilibrium dynamics appear to be independent of initial state variables
species_pair.parms <- data.frame(r1 = r_zoops_estimate, 
                                r2 = r_benthos_estimate, 
                                K1 = K_zoops_estimate, 
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

species_pair_dynamics <- safe.runsteady(y = stickleback.init.states, func = ECD_model, parms = stickleback.parms, stol = 1e-5)
species_pair_dynamics$y # equilibrium densities

# estimates at first glance seem reasonable. Biomass of consumers is less than their prey, which is what we would
# expect based on trophic transfer efficiency
matplot.0D(ode(species_pair.init.states, 1:1000, ECD_model, species_pair.parms))


## Solitary lakes
solitary.init.states <- c(R1=1, R2=1, C1=0.5) # note that equilibrium dynamics appear to be independent of initial state variables
solitary.parms <- data.frame(r1 = r_zoops_estimate, 
                             r2 = r_benthos_estimate, 
                             K1 = K_zoops_estimate, 
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
                             m1 = m_solitary_estimate) 
                             #m2 = m_benthic_estimate)

solitary_dynamics <- safe.runsteady(y = solitary.init.states, func = ECD_model.C1.3sp, parms = solitary.parms, stol = 1e-5)
solitary_dynamics$y # equilibrium densities

# estimates at first glance seem reasonable. Biomass of consumers is less than their prey, which is what we would
# expect based on trophic transfer efficiency
matplot.0D(ode(solitary.init.states, 1:1000, ECD_model.C1.3sp, solitary.parms))


## compare species-pair and solitary lakes
solitary_total_R <- solitary_dynamics$y["R1"] + solitary_dynamics$y["R2"]; solitary_total_R
species_pair_total_R <- species_pair_dynamics$y["R1"] + species_pair_dynamics$y["R2"]; species_pair_total_R 

species_pair_total_R/solitary_total_R - 1 # total resource abundances are lower in species pair lakes

species_pair_dynamics$y["R1"]/solitary_dynamics$y["R1"] - 1 # zooplankton densities are lower in species pair lakes
