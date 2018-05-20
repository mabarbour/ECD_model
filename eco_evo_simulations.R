
#### PRELIMINARIES ####

## Source in required functions
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')
source('eco_evo_sim_function.R')
source('Consumer_Resource_Functions.R')
source('CR_dynamics.R')

#### GENERAL PARAMETERS AND STATE VARIABLES #### 

## Intrinsic growth rates
r <- 1
r1 <- r 
r2 <- r 

## Carrying capacity
K <- 4
K1 <- K
K2 <- K

## Conversion efficiency
e <- 0.8
e11 <- e 
e12 <- e 
e21 <- e 
e22 <- e 

## Mortality rate
m <- 1
m1 <- m
m2 <- m

## General, non-evolving parameter set for 2C,2R model
general_parameters <- data.frame(r1 = r1, r2 = r2, K1 = K1, K2 = K2, e11 = e11, e12 = e12, e21 = e21, e22 = e22, m1 = m1, m2 = m2)

## Attack rate
A <- 2
aii <- A*0.6

## Habitat preference (for LawlorSmith and McCann models)
w <- 0.6
wii <- w
wij <- (1 - w)

## Handling time (only for McCann model)
h <- 0.4
h11 <- h
h12 <- h
h21 <- h
h22 <- h

## State variables
R <- 2
R1 <- R
R2 <- R
C <- 1
C1 <- C
C2 <- C

states_2C_2R <- data.frame(R1 = R1, R2 = R2, C1 = C1, C2 = C2)
states_1C_2R <- data.frame(R1 = R1, R2 = R2, C1 = C1)

## Trade-offs. Borrowed from Sargent and Otto, 2006 Am Nat
a12_tradeoff <- expression(A*(1-(a11m/A)^n)^(1/n)) 
a21_tradeoff <- expression(A*(1-(a22m/A)^n)^(1/n)) 

n_concaveDN <- 1.15
n_linear <- 1
n_concaveUP <- 0.85

concaveDN_df <- data.frame(n = n_concaveDN, aii = seq(0,A,0.01)) %>% mutate(aij = A*(1-(aii/A)^n_concaveDN)^(1/n_concaveDN))
linear_df <- data.frame(n = n_linear, aii = seq(0,A,0.01)) %>% mutate(aij = A*(1-(aii/A)^n_linear)^(1/n_linear))
concaveUP_df <- data.frame(n = n_concaveUP, aii = seq(0,A,0.01)) %>% mutate(aij = A*(1-(aii/A)^n_concaveUP)^(1/n_concaveUP))

bind_rows(concaveDN_df, linear_df, concaveUP_df) %>%
  ggplot(., aes(x=aii, y=aij, color = factor(n))) + geom_line() + ggtitle("Attack rate trade-offs")

bind_rows(concaveDN_df, linear_df, concaveUP_df) %>%
  ggplot(., aes(x=aii/(aii+aij), y=aii+aij, color = factor(n))) + geom_line() + ggtitle("Relationship between specialization and total attack rate")

######################################### MACARTHUR MODEL #######################################################

## Mutant growth rates
mC1_MacArthur <- expression(e11*a11m*R1 + e12*a12m*R2 - m1)
mC2_MacArthur <- expression(e21*a21m*R1 + e22*a22m*R2 - m2)

#### PARAMETERS ####

## Linear
params_2C_2R_MacArthur_linear <- data.frame(general_parameters,
                                            a11 = aii, 
                                            a12 = eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_linear)),                     
                                            a21 = eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_linear)), 
                                            a22 = aii)                      
params_1C_2R_MacArthur_linear <- params_2C_2R_MacArthur_linear[c("r1","r2","K1","K2","e11","e12","a11","a12","m1")]


## concaveDN
params_2C_2R_MacArthur_concaveDN <- params_2C_2R_MacArthur_linear
params_2C_2R_MacArthur_concaveDN["a12"] <- eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_concaveDN))
params_2C_2R_MacArthur_concaveDN["a21"] <- eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_concaveDN))
params_1C_2R_MacArthur_concaveDN <- params_2C_2R_MacArthur_concaveDN[c("r1","r2","K1","K2","e11","e12","a11","a12","m1")]

## concaveUP
params_2C_2R_MacArthur_concaveUP <- params_2C_2R_MacArthur_linear
params_2C_2R_MacArthur_concaveUP["a12"] <- eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_concaveUP))
params_2C_2R_MacArthur_concaveUP["a21"] <- eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_concaveUP))
params_1C_2R_MacArthur_concaveUP <- params_2C_2R_MacArthur_concaveUP[c("r1","r2","K1","K2","e11","e12","a11","a12","m1")]

#### SIMULATIONS ####

## Linear
sim_2C_2R_MacArthur_linear <- eco_evo_CR(init_parameters = params_2C_2R_MacArthur_linear, 
                                init_states = states_2C_2R, 
                                eco_CR_model = MacArthur_2C_2R, 
                                AD_CR_models = c(mC1_MacArthur, mC2_MacArthur),
                                species_traits = c(1,2), 
                                eco_param_names = c("a11","a22"), 
                                eco_tradeoff_names = c("a12","a21"),  
                                evo_param_names = c("a11m","a22m"),  
                                evo_tradeoff_names = c("a12m","a21m"), 
                                tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                extra_tradeoff_params = c(A=A, n=n_linear))
sim_1C_2R_MacArthur_linear <- eco_evo_CR(init_parameters = params_1C_2R_MacArthur_linear, 
                                init_states = states_1C_2R, 
                                eco_CR_model = MacArthur_1C_2R, 
                                AD_CR_models = c(mC1_MacArthur),
                                species_traits = c(1), 
                                eco_param_names = c("a11"), 
                                eco_tradeoff_names = c("a12"),  
                                evo_param_names = c("a11m"),  
                                evo_tradeoff_names = c("a12m"), 
                                tradeoff_functions = c(a12_tradeoff),
                                extra_tradeoff_params = c(A=A, n=n_linear))

## concaveDN
sim_2C_2R_MacArthur_concaveDN <- eco_evo_CR(init_parameters = params_2C_2R_MacArthur_concaveDN, 
                                init_states = states_2C_2R, 
                                eco_CR_model = MacArthur_2C_2R, 
                                AD_CR_models = c(mC1_MacArthur, mC2_MacArthur),
                                species_traits = c(1,2), 
                                eco_param_names = c("a11","a22"), 
                                eco_tradeoff_names = c("a12","a21"),  
                                evo_param_names = c("a11m","a22m"),  
                                evo_tradeoff_names = c("a12m","a21m"), 
                                tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                extra_tradeoff_params = c(A=A, n=n_concaveDN))
sim_1C_2R_MacArthur_concaveDN <- eco_evo_CR(init_parameters = params_1C_2R_MacArthur_concaveDN, 
                                init_states = states_1C_2R, 
                                eco_CR_model = MacArthur_1C_2R, 
                                AD_CR_models = c(mC1_MacArthur),
                                species_traits = c(1), 
                                eco_param_names = c("a11"), 
                                eco_tradeoff_names = c("a12"),  
                                evo_param_names = c("a11m"),  
                                evo_tradeoff_names = c("a12m"), 
                                tradeoff_functions = c(a12_tradeoff),
                                extra_tradeoff_params = c(A=A, n=n_concaveDN))
## concaveUP
sim_2C_2R_MacArthur_concaveUP <- eco_evo_CR(init_parameters = params_2C_2R_MacArthur_concaveUP, 
                                       init_states = states_2C_2R, 
                                       eco_CR_model = MacArthur_2C_2R, 
                                       AD_CR_models = c(mC1_MacArthur, mC2_MacArthur),
                                       species_traits = c(1,2), 
                                       eco_param_names = c("a11","a22"), 
                                       eco_tradeoff_names = c("a12","a21"),  
                                       evo_param_names = c("a11m","a22m"),  
                                       evo_tradeoff_names = c("a12m","a21m"), 
                                       tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                       extra_tradeoff_params = c(A=A, n=n_concaveUP))
sim_1C_2R_MacArthur_concaveUP <- eco_evo_CR(init_parameters = params_1C_2R_MacArthur_concaveUP, 
                                       init_states = states_1C_2R, 
                                       eco_CR_model = MacArthur_1C_2R, 
                                       AD_CR_models = c(mC1_MacArthur),
                                       species_traits = c(1), 
                                       eco_param_names = c("a11"), 
                                       eco_tradeoff_names = c("a12"),  
                                       evo_param_names = c("a11m"),  
                                       evo_tradeoff_names = c("a12m"), 
                                       tradeoff_functions = c(a12_tradeoff),
                                       extra_tradeoff_params = c(A=A, n=n_concaveUP))

######################################### LAWLOR SMITH COARSE-GRAIN MODEL #######################################################


## Mutant growth rates
mC1_LawlorSmith <- expression(e11*w11*a11m*R1 + e12*w12*a12m*R2 - m1)
mC2_LawlorSmith <- expression(e21*w21*a21m*R1 + e22*w22*a22m*R2 - m2)

#### PARAMETERS ####

## Linear
params_2C_2R_LawlorSmith_linear <- data.frame(general_parameters,
                                              a11 = aii,
                                              a12 = eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_linear)), 
                                              a21 = eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_linear)), 
                                              a22 = aii,
                                              w11 = wii, 
                                              w12 = wij, 
                                              w21 = wij,
                                              w22 = wii)
params_1C_2R_LawlorSmith_linear <- params_2C_2R_LawlorSmith_linear[c("r1","r2","K1","K2","e11","e12","a11","a12","m1","w11","w12","w21","w22")]


## concaveDN
params_2C_2R_LawlorSmith_concaveDN <- params_2C_2R_LawlorSmith_linear
params_2C_2R_LawlorSmith_concaveDN["a12"] <- eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_concaveDN))
params_2C_2R_LawlorSmith_concaveDN["a21"] <- eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_concaveDN))
params_1C_2R_LawlorSmith_concaveDN <- params_2C_2R_LawlorSmith_concaveDN[c("r1","r2","K1","K2","e11","e12","a11","a12","m1","w11","w12","w21","w22")]

## concaveUP
params_2C_2R_LawlorSmith_concaveUP <- params_2C_2R_LawlorSmith_linear
params_2C_2R_LawlorSmith_concaveUP["a12"] <- eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_concaveUP))
params_2C_2R_LawlorSmith_concaveUP["a21"] <- eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_concaveUP))
params_1C_2R_LawlorSmith_concaveUP <- params_2C_2R_LawlorSmith_concaveUP[c("r1","r2","K1","K2","e11","e12","a11","a12","m1","w11","w12","w21","w22")]

#### SIMULATIONS ####

## Linear
sim_2C_2R_LawlorSmith_linear <- eco_evo_CR(init_parameters = params_2C_2R_LawlorSmith_linear, 
                                       init_states = states_2C_2R, 
                                       eco_CR_model = LawlorSmith_2C_2R, 
                                       AD_CR_models = c(mC1_LawlorSmith, mC2_LawlorSmith),
                                       species_traits = c(1,2), 
                                       eco_param_names = c("a11","a22"), 
                                       eco_tradeoff_names = c("a12","a21"),  
                                       evo_param_names = c("a11m","a22m"),  
                                       evo_tradeoff_names = c("a12m","a21m"), 
                                       tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                       extra_tradeoff_params = c(A=A, n=n_linear))
sim_1C_2R_LawlorSmith_linear <- eco_evo_CR(init_parameters = params_1C_2R_LawlorSmith_linear, 
                                       init_states = states_1C_2R, 
                                       eco_CR_model = LawlorSmith_1C_2R, 
                                       AD_CR_models = c(mC1_LawlorSmith),
                                       species_traits = c(1), 
                                       eco_param_names = c("a11"), 
                                       eco_tradeoff_names = c("a12"),  
                                       evo_param_names = c("a11m"),  
                                       evo_tradeoff_names = c("a12m"), 
                                       tradeoff_functions = c(a12_tradeoff),
                                       extra_tradeoff_params = c(A=A, n=n_linear))

## concaveDN
# max mut iter reach!
sim_2C_2R_LawlorSmith_concaveDN <- eco_evo_CR(init_parameters = params_2C_2R_LawlorSmith_concaveDN, 
                                       init_states = states_2C_2R, 
                                       eco_CR_model = LawlorSmith_2C_2R, 
                                       AD_CR_models = c(mC1_LawlorSmith, mC2_LawlorSmith),
                                       species_traits = c(1,2), 
                                       eco_param_names = c("a11","a22"), 
                                       eco_tradeoff_names = c("a12","a21"),  
                                       evo_param_names = c("a11m","a22m"),  
                                       evo_tradeoff_names = c("a12m","a21m"), 
                                       tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                       extra_tradeoff_params = c(A=A, n=n_concaveDN))
# note that although the max mutation iterations are reached, the evolving parameters appear to oscillate near high values, i.e. they are effectively at an ESS
plot(sim_2C_2R_LawlorSmith_concaveDN[,c("a11","a22")]) 
plot(sim_2C_2R_LawlorSmith_concaveDN[,c("a11","a22")], xlim=c(1.85,1.95)) 

sim_1C_2R_LawlorSmith_concaveDN <- eco_evo_CR(init_parameters = params_1C_2R_LawlorSmith_concaveDN, 
                                       init_states = states_1C_2R, 
                                       eco_CR_model = LawlorSmith_1C_2R, 
                                       AD_CR_models = c(mC1_LawlorSmith),
                                       species_traits = c(1), 
                                       eco_param_names = c("a11"), 
                                       eco_tradeoff_names = c("a12"),  
                                       evo_param_names = c("a11m"),  
                                       evo_tradeoff_names = c("a12m"), 
                                       tradeoff_functions = c(a12_tradeoff),
                                       extra_tradeoff_params = c(A=A, n=n_concaveDN))
## concaveUP
sim_2C_2R_LawlorSmith_concaveUP <- eco_evo_CR(init_parameters = params_2C_2R_LawlorSmith_concaveUP, 
                                        init_states = states_2C_2R, 
                                        eco_CR_model = LawlorSmith_2C_2R, 
                                        AD_CR_models = c(mC1_LawlorSmith, mC2_LawlorSmith),
                                        species_traits = c(1,2), 
                                        eco_param_names = c("a11","a22"), 
                                        eco_tradeoff_names = c("a12","a21"),  
                                        evo_param_names = c("a11m","a22m"),  
                                        evo_tradeoff_names = c("a12m","a21m"), 
                                        tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                        extra_tradeoff_params = c(A=A, n=n_concaveUP))
sim_1C_2R_LawlorSmith_concaveUP <- eco_evo_CR(init_parameters = params_1C_2R_LawlorSmith_concaveUP, 
                                        init_states = states_1C_2R, 
                                        eco_CR_model = LawlorSmith_1C_2R, 
                                        AD_CR_models = c(mC1_LawlorSmith),
                                        species_traits = c(1), 
                                        eco_param_names = c("a11"), 
                                        eco_tradeoff_names = c("a12"),  
                                        evo_param_names = c("a11m"),  
                                        evo_tradeoff_names = c("a12m"), 
                                        tradeoff_functions = c(a12_tradeoff),
                                        extra_tradeoff_params = c(A=A, n=n_concaveUP))

######################################### MCCANN MODEL #######################################################


## Mutant growth rates
mC1_McCann <- expression(e11*(((w11 * R1)/(w11 * R1 + w12 * R2)) * a11m * R1)/(1 + ((w11 * R1)/(w11 * R1 + w12 * R2)) * a11m * h11 * R1 + ((w12 * R2)/(w11 * R1 + w12 * R2)) * a12m * h12 * R2) + e12*(((w12 * R2)/(w11 * R1 + w12 * R2)) * a12m * R2)/(1 + ((w11 * R1)/(w11 * R1 + w12 * R2)) * a11m * h11 * R1 + ((w12 * R2)/(w11 * R1 + w12 * R2)) * a12m * h12 * R2) - m1)
mC2_McCann <- expression(e21*(((w21 * R1)/(w21 * R1 + w22 * R2)) * a21m * R1)/(1 + ((w21 * R1)/(w21 * R1 + w22 * R2)) * a21m * h21 * R1 + ((w22 * R2)/(w21 * R1 + w22 * R2)) * a22m * h22 * R2) + e22*(((w22 * R2)/(w21 * R1 + w22 * R2)) * a22m * R2)/(1 + ((w21 * R1)/(w21 * R1 + w22 * R2)) * a21m * h21 * R1 + ((w22 * R2)/(w21 * R1 + w22 * R2)) * a22m * h22 * R2) - m2)


#### PARAMETERS ####

## Linear
params_2C_2R_McCann_linear <- data.frame(general_parameters,
                                         a11 = aii,
                                         a12 = eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_linear)), 
                                         a21 = eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_linear)), 
                                         a22 = aii,
                                         w11 = wii, 
                                         w12 = wij, 
                                         w21 = wij,
                                         w22 = wii,
                                         h11 = h11,
                                         h12 = h12,
                                         h21 = h21,
                                         h22 = h22)
                                            
params_1C_2R_McCann_linear <- params_2C_2R_McCann_linear[c("r1","r2","K1","K2","e11","e12","a11","a12","m1","w11","w12","w21","w22","h11","h12","h21","h22")]


## concaveDN
params_2C_2R_McCann_concaveDN <- params_2C_2R_McCann_linear
params_2C_2R_McCann_concaveDN["a12"] <- eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_concaveDN))
params_2C_2R_McCann_concaveDN["a21"] <- eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_concaveDN))
params_1C_2R_McCann_concaveDN <- params_2C_2R_McCann_concaveDN[c("r1","r2","K1","K2","e11","e12","a11","a12","m1","w11","w12","w21","w22","h11","h12","h21","h22")]

## concaveUP
params_2C_2R_McCann_concaveUP <- params_2C_2R_McCann_linear
params_2C_2R_McCann_concaveUP["a12"] <- eval(a12_tradeoff, data.frame(A=A, a11m=aii, n=n_concaveUP))
params_2C_2R_McCann_concaveUP["a21"] <- eval(a21_tradeoff, data.frame(A=A, a22m=aii, n=n_concaveUP))
params_1C_2R_McCann_concaveUP <- params_2C_2R_McCann_concaveUP[c("r1","r2","K1","K2","e11","e12","a11","a12","m1","w11","w12","w21","w22","h11","h12","h21","h22")]

#### SIMULATIONS ####

## Linear
sim_2C_2R_McCann_linear <- eco_evo_CR(init_parameters = params_2C_2R_McCann_linear, 
                                         init_states = states_2C_2R, 
                                         eco_CR_model = McCann_2C_2R, 
                                         AD_CR_models = c(mC1_McCann, mC2_McCann),
                                         species_traits = c(1,2), 
                                         eco_param_names = c("a11","a22"), 
                                         eco_tradeoff_names = c("a12","a21"),  
                                         evo_param_names = c("a11m","a22m"),  
                                         evo_tradeoff_names = c("a12m","a21m"), 
                                         tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                         extra_tradeoff_params = c(A=A, n=n_linear))
sim_1C_2R_McCann_linear <- eco_evo_CR(init_parameters = params_1C_2R_McCann_linear, 
                                         init_states = states_1C_2R, 
                                         eco_CR_model = McCann_1C_2R, 
                                         AD_CR_models = c(mC1_McCann),
                                         species_traits = c(1), 
                                         eco_param_names = c("a11"), 
                                         eco_tradeoff_names = c("a12"),  
                                         evo_param_names = c("a11m"),  
                                         evo_tradeoff_names = c("a12m"), 
                                         tradeoff_functions = c(a12_tradeoff),
                                         extra_tradeoff_params = c(A=A, n=n_linear))

## concaveDN
# max mut iter
sim_2C_2R_McCann_concaveDN <- eco_evo_CR(init_parameters = params_2C_2R_McCann_concaveDN, 
                                         init_states = states_2C_2R, 
                                         eco_CR_model = McCann_2C_2R, 
                                         AD_CR_models = c(mC1_McCann, mC2_McCann),
                                         species_traits = c(1,2), 
                                         eco_param_names = c("a11","a22"), 
                                         eco_tradeoff_names = c("a12","a21"),  
                                         evo_param_names = c("a11m","a22m"),  
                                         evo_tradeoff_names = c("a12m","a21m"), 
                                         tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                         extra_tradeoff_params = c(A=A, n=n_concaveDN))
# note that although the max mutation iterations are reached, the evolving parameters appear to oscillate near high values. i.e. they are effectively at an ESS
plot(sim_2C_2R_McCann_concaveDN[,c("a11","a22")]) 
plot(sim_2C_2R_McCann_concaveDN[,c("a11","a22")], xlim=c(1.85,1.95)) 

sim_1C_2R_McCann_concaveDN <- eco_evo_CR(init_parameters = params_1C_2R_McCann_concaveDN, 
                                         init_states = states_1C_2R, 
                                         eco_CR_model = McCann_1C_2R, 
                                         AD_CR_models = c(mC1_McCann),
                                         species_traits = c(1), 
                                         eco_param_names = c("a11"), 
                                         eco_tradeoff_names = c("a12"),  
                                         evo_param_names = c("a11m"),  
                                         evo_tradeoff_names = c("a12m"), 
                                         tradeoff_functions = c(a12_tradeoff),
                                         extra_tradeoff_params = c(A=A, n=n_concaveDN))
## concaveUP
sim_2C_2R_McCann_concaveUP <- eco_evo_CR(init_parameters = params_2C_2R_McCann_concaveUP, 
                                          init_states = states_2C_2R, 
                                          eco_CR_model = McCann_2C_2R, 
                                          AD_CR_models = c(mC1_McCann, mC2_McCann),
                                          species_traits = c(1,2), 
                                          eco_param_names = c("a11","a22"), 
                                          eco_tradeoff_names = c("a12","a21"),  
                                          evo_param_names = c("a11m","a22m"),  
                                          evo_tradeoff_names = c("a12m","a21m"), 
                                          tradeoff_functions = c(a12_tradeoff, a21_tradeoff),
                                          extra_tradeoff_params = c(A=A, n=n_concaveUP))
sim_1C_2R_McCann_concaveUP <- eco_evo_CR(init_parameters = params_1C_2R_McCann_concaveUP, 
                                          init_states = states_1C_2R, 
                                          eco_CR_model = McCann_1C_2R, 
                                          AD_CR_models = c(mC1_McCann),
                                          species_traits = c(1), 
                                          eco_param_names = c("a11"), 
                                          eco_tradeoff_names = c("a12"),  
                                          evo_param_names = c("a11m"),  
                                          evo_tradeoff_names = c("a12m"), 
                                          tradeoff_functions = c(a12_tradeoff),
                                          extra_tradeoff_params = c(A=A, n=n_concaveUP))


#### PLOTS ####

## Gather data. Note different formula for effective attack rate of consumer 1 "a_eff_C1". This accounts for the differences between fine- and coarse-grained environments
sim_data <- bind_rows(
  mutate(sim_2C_2R_MacArthur_linear, Model = "MacArthur", n = n_linear, Competitor = "Yes", a_eff_C1 = a11+a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_MacArthur_linear, Model = "MacArthur", n = n_linear, Competitor = "No", a_eff_C1 = a11+a12, Ctot = C1),
  mutate(sim_2C_2R_MacArthur_concaveDN, Model = "MacArthur", n = n_concaveDN, Competitor = "Yes", a_eff_C1 = a11+a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_MacArthur_concaveDN, Model = "MacArthur", n = n_concaveDN, Competitor = "No", a_eff_C1 = a11+a12, Ctot = C1),
  mutate(sim_2C_2R_MacArthur_concaveUP, Model = "MacArthur", n = n_concaveUP, Competitor = "Yes", a_eff_C1 = a11+a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_MacArthur_concaveUP, Model = "MacArthur", n = n_concaveUP, Competitor = "No", a_eff_C1 = a11+a12, Ctot = C1),
  mutate(sim_2C_2R_LawlorSmith_linear, Model = "LawlorSmith", n = n_linear, Competitor = "Yes", a_eff_C1 = w11*a11+w12*a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_LawlorSmith_linear, Model = "LawlorSmith", n = n_linear, Competitor = "No", a_eff_C1 = w11*a11+w12*a12, Ctot = C1),
  mutate(sim_2C_2R_LawlorSmith_concaveDN, Model = "LawlorSmith", n = n_concaveDN, Competitor = "Yes", a_eff_C1 = w11*a11+w12*a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_LawlorSmith_concaveDN, Model = "LawlorSmith", n = n_concaveDN, Competitor = "No", a_eff_C1 = w11*a11+w12*a12, Ctot = C1),
  mutate(sim_2C_2R_LawlorSmith_concaveUP, Model = "LawlorSmith", n = n_concaveUP, Competitor = "Yes", a_eff_C1 = w11*a11+w12*a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_LawlorSmith_concaveUP, Model = "LawlorSmith", n = n_concaveUP, Competitor = "No", a_eff_C1 = w11*a11+w12*a12, Ctot = C1),
  mutate(sim_2C_2R_McCann_linear, Model = "McCann", n = n_linear, Competitor = "Yes", a_eff_C1 = w11*a11+w12*a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_McCann_linear, Model = "McCann", n = n_linear, Competitor = "No", a_eff_C1 = w11*a11+w12*a12, Ctot = C1),
  mutate(sim_2C_2R_McCann_concaveDN, Model = "McCann", n = n_concaveDN, Competitor = "Yes", a_eff_C1 = w11*a11+w12*a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_McCann_concaveDN, Model = "McCann", n = n_concaveDN, Competitor = "No", a_eff_C1 = w11*a11+w12*a12, Ctot = C1),
  mutate(sim_2C_2R_McCann_concaveUP, Model = "McCann", n = n_concaveUP, Competitor = "Yes", a_eff_C1 = w11*a11+w12*a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_McCann_concaveUP, Model = "McCann", n = n_concaveUP, Competitor = "No", a_eff_C1 = w11*a11+w12*a12, Ctot = C1)) 


## Trait Divergence
#ggplot(sim_data, aes(x=sequence, y=a11/(a11+a12), color=Trade_off, linetype=Competitor)) + 
#  geom_line() + facet_wrap(~Model) + geom_hline(yintercept=0.5, linetype="dotted") + ggtitle("Trait Divergence")

## Effective Attack Rate
#ggplot(sim_data, aes(x=sequence, y=a_eff_C1, color=Trade_off, linetype=Competitor)) + 
#  geom_line() + facet_wrap(~Model) + geom_hline(yintercept=0.5, linetype="dotted") + ggtitle("Effective attack rate")

## Abundances
#gather_abundance <- sim_data %>%
#  select(sequence, Model, Trade_off, Competitor, R1, R2, C1, C2) %>%
#  gather(Species, Density, R1:C2)

#ggplot(gather_abundance %>% filter(Species %in% "R1"), aes(x=sequence, y=Density, color=Species, linetype=Competitor)) +
#  geom_line() + facet_grid(Trade_off~Model)

#ggplot(gather_abundance %>% filter(Species %in% "R2"), aes(x=sequence, y=Density, color=Species, linetype=Competitor)) +
#  geom_line() + facet_grid(Trade_off~Model)


## Stability
#ggplot(sim_data, aes(x=sequence, y=-1*max.Re.eigen, color=Trade_off, linetype=Competitor)) + 
#  geom_line() + facet_wrap(~Model) + ggtitle("Stability")


#### PREDICTING EFFECTS OF CHARACTER DISPLACEMENT ON C-R DYNAMICS ####

tradeoff_df <- bind_rows(mutate(concaveDN_df, n = n_concaveDN), 
                         mutate(linear_df, n = n_linear), 
                         mutate(concaveUP_df, n = n_concaveUP)) %>%
  data.frame(., general_parameters) %>%
  mutate(a11 = aii, a12 = aij, a21 = aij, a22 = aii,
         w11 = wii, w12 = wij, w21 = wij, w22 = wii)

## MacArthur
MacArthur_displacements <- bind_rows(
  mutate(sim_2C_2R_MacArthur_linear, Model = "MacArthur", n = n_linear, Competitor = "Yes")[c(1,nrow(sim_2C_2R_MacArthur_linear)), ], 
  mutate(sim_1C_2R_MacArthur_linear, Model = "MacArthur", n = n_linear, Competitor = "No")[c(1,nrow(sim_1C_2R_MacArthur_linear)), ],
  mutate(sim_2C_2R_MacArthur_concaveDN, Model = "MacArthur", n = n_concaveDN, Competitor = "Yes")[c(1,nrow(sim_2C_2R_MacArthur_concaveDN)), ], 
  mutate(sim_1C_2R_MacArthur_concaveDN, Model = "MacArthur", n = n_concaveDN, Competitor = "No")[c(1,nrow(sim_1C_2R_MacArthur_concaveDN)), ],
  mutate(sim_2C_2R_MacArthur_concaveUP, Model = "MacArthur", n = n_concaveUP, Competitor = "Yes")[c(1,nrow(sim_2C_2R_MacArthur_concaveUP)), ], 
  mutate(sim_1C_2R_MacArthur_concaveUP, Model = "MacArthur", n = n_concaveUP, Competitor = "No")[c(1,nrow(sim_1C_2R_MacArthur_concaveUP)), ]) 
MacArthur_displacements$Time <- rep(c("Begin","End"),6)

ggplot(tradeoff_df, aes(x=aii/(aii+aij), y=aii+aij, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=a11+a12, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

MacArthur_2C_2R_tradeoff_dynamics <- CR_dynamics(init_parameters = select(tradeoff_df, r1:a22), init_states = states_2C_2R, eco_CR_model = MacArthur_2C_2R)
MacArthur_1C_2R_tradeoff_dynamics <- CR_dynamics(init_parameters = select(tradeoff_df, r1:a12), init_states = states_1C_2R, eco_CR_model = MacArthur_1C_2R)

# whoa...I can predict the dynamics of both the 3 species and 4 species system with just the abundance data from the 4 species system...This must be because when only 1C is present, then the consumer goes back to complete generalist...
cbind.data.frame(tradeoff_df, MacArthur_2C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=R1+R2, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=R1+R2, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

cbind.data.frame(tradeoff_df, MacArthur_1C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=R1+R2, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=R1+R2, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

cbind.data.frame(tradeoff_df, MacArthur_2C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=-1*max.Re.eigen, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=-1*max.Re.eigen, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

cbind.data.frame(tradeoff_df, MacArthur_1C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=-1*max.Re.eigen, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=-1*max.Re.eigen, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

## Lawlor Smith
LawlorSmith_displacements <- bind_rows(
  mutate(sim_2C_2R_LawlorSmith_linear, Model = "LawlorSmith", n = n_linear, Competitor = "Yes")[c(1,nrow(sim_2C_2R_LawlorSmith_linear)), ], 
  mutate(sim_1C_2R_LawlorSmith_linear, Model = "LawlorSmith", n = n_linear, Competitor = "No")[c(1,nrow(sim_1C_2R_LawlorSmith_linear)), ],
  mutate(sim_2C_2R_LawlorSmith_concaveDN, Model = "LawlorSmith", n = n_concaveDN, Competitor = "Yes")[c(1,nrow(sim_2C_2R_LawlorSmith_concaveDN)), ], 
  mutate(sim_1C_2R_LawlorSmith_concaveDN, Model = "LawlorSmith", n = n_concaveDN, Competitor = "No")[c(1,nrow(sim_1C_2R_LawlorSmith_concaveDN)), ],
  mutate(sim_2C_2R_LawlorSmith_concaveUP, Model = "LawlorSmith", n = n_concaveUP, Competitor = "Yes")[c(1,nrow(sim_2C_2R_LawlorSmith_concaveUP)), ], 
  mutate(sim_1C_2R_LawlorSmith_concaveUP, Model = "LawlorSmith", n = n_concaveUP, Competitor = "No")[c(1,nrow(sim_1C_2R_LawlorSmith_concaveUP)), ]) 
LawlorSmith_displacements$Time <- rep(c("Begin","End"),6)

ggplot(tradeoff_df, aes(x=aii/(aii+aij), y=wii*aii+wij*aij, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(LawlorSmith_displacements, Time=="End"), aes(x = a11/(a11+a12), y=w11*a11+w12*a12, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

LawlorSmith_2C_2R_tradeoff_dynamics <- CR_dynamics(init_parameters = select(tradeoff_df, r1:w22), init_states = states_2C_2R, eco_CR_model = LawlorSmith_2C_2R)
LawlorSmith_1C_2R_tradeoff_dynamics <- CR_dynamics(init_parameters = select(tradeoff_df, r1:a12, w11, w12), init_states = states_1C_2R, eco_CR_model = LawlorSmith_1C_2R)

cbind.data.frame(tradeoff_df, LawlorSmith_2C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=R1+R2, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(LawlorSmith_displacements, Time=="End"), aes(x = a11/(a11+a12), y=R1+R2, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

cbind.data.frame(tradeoff_df, LawlorSmith_1C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=R1+R2, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(LawlorSmith_displacements, Time=="End"), aes(x = a11/(a11+a12), y=R1+R2, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

cbind.data.frame(tradeoff_df, LawlorSmith_2C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=-1*max.Re.eigen, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(LawlorSmith_displacements, Time=="End"), aes(x = a11/(a11+a12), y=-1*max.Re.eigen, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

cbind.data.frame(tradeoff_df, LawlorSmith_1C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=-1*max.Re.eigen, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(LawlorSmith_displacements, Time=="End"), aes(x = a11/(a11+a12), y=-1*max.Re.eigen, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))


## McCann
McCann_displacements <- bind_rows(
  mutate(sim_2C_2R_McCann_linear, Model = "McCann", n = n_linear, Competitor = "Yes")[c(1,nrow(sim_2C_2R_McCann_linear)), ], 
  mutate(sim_1C_2R_McCann_linear, Model = "McCann", n = n_linear, Competitor = "No")[c(1,nrow(sim_1C_2R_McCann_linear)), ],
  mutate(sim_2C_2R_McCann_concaveDN, Model = "McCann", n = n_concaveDN, Competitor = "Yes")[c(1,nrow(sim_2C_2R_McCann_concaveDN)), ], 
  mutate(sim_1C_2R_McCann_concaveDN, Model = "McCann", n = n_concaveDN, Competitor = "No")[c(1,nrow(sim_1C_2R_McCann_concaveDN)), ],
  mutate(sim_2C_2R_McCann_concaveUP, Model = "McCann", n = n_concaveUP, Competitor = "Yes")[c(1,nrow(sim_2C_2R_McCann_concaveUP)), ], 
  mutate(sim_1C_2R_McCann_concaveUP, Model = "McCann", n = n_concaveUP, Competitor = "No")[c(1,nrow(sim_1C_2R_McCann_concaveUP)), ]) 
McCann_displacements$Time <- rep(c("Begin","End"),6)

ggplot(tradeoff_df, aes(x=aii/(aii+aij), y=wii*aii+wij*aij, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(McCann_displacements, Time=="End"), aes(x = a11/(a11+a12), y=w11*a11+w12*a12, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

McCann_2C_2R_tradeoff_dynamics <- CR_dynamics(init_parameters = select(tradeoff_df, r1:w22), init_states = states_2C_2R, eco_CR_model = McCann_2C_2R)
McCann_1C_2R_tradeoff_dynamics <- CR_dynamics(init_parameters = select(tradeoff_df, r1:a12, w11, w12), init_states = states_1C_2R, eco_CR_model = McCann_1C_2R)

cbind.data.frame(tradeoff_df, McCann_2C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=R1+R2, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(McCann_displacements, Time=="End"), aes(x = a11/(a11+a12), y=R1+R2, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

cbind.data.frame(tradeoff_df, McCann_1C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=R1+R2, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(McCann_displacements, Time=="End"), aes(x = a11/(a11+a12), y=R1+R2, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

# appears that in some case, stability would still increase with displacement in the McCann model...(why?)
cbind.data.frame(tradeoff_df, McCann_2C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=-1*max.Re.eigen, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(McCann_displacements, Time=="End"), aes(x = a11/(a11+a12), y=-1*max.Re.eigen, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))

cbind.data.frame(tradeoff_df, McCann_1C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=aii/(aii+aij), y=-1*max.Re.eigen, linetype = factor(n))) + 
  geom_line() +
  geom_point(data = filter(McCann_displacements, Time=="End"), aes(x = a11/(a11+a12), y=-1*max.Re.eigen, fill=Competitor), shape=21, size=10, inherit.aes = FALSE) +
  scale_fill_manual(values=c("white","black"))


