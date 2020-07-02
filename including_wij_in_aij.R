## Response to reviewer #2

# In this script, I address reviewer #2's question:
# How would Fig. 3(d,e,f) change if eq. (2) was used with 
# tilde(a)_ij=w_ij a_ij (and the trade-off applied to tilde(a)_ij?).


## Load required libraries
library(deSolve) # numerical integration library
library(rootSolve) # for runsteady integration function
library(pse) # for latin hypercube
library(plyr) # for ldply
library(dplyr) # for manipulating data
library(tidyr) # for tidying data
library(cowplot) # for better base, ggplot graphics

## color palette
cbPalette <- c("#D55E00", "#0072B2", "#009E73")

## set plot theme
theme_set(theme_cowplot())

#### GET REQUIRED FUNCTIONS ####

source('identify_steady_state.R')
source('eco_evo_sim_function.R')
source('Consumer_Resource_Functions.R')
source('CR_dynamics.R')

#### GENERAL PARAMETERS AND STATE VARIABLES #### 

## Intrinsic growth rates
r <- 1
r1 <- r + 0.1
r2 <- r 

## Carrying capacity
K <- 4
K1 <- K - 0.1
K2 <- K

## Conversion efficiency
e <- 0.8
e11 <- e + 0.1
e12 <- e - 0.1
e21 <- e 
e22 <- e 

## Mortality rate
m <- 1
m1 <- m
m2 <- m

## Habitat preference (for LawlorSmith and McCann models)
w <- 0.61 # adjusted to 0.61 from 0.6 to maintain consistency with previous figure
wii <- w
wij <- (1 - w)

## General, non-evolving parameter set for 2C,2R model
general_parameters <- data.frame(r1 = r1, r2 = r2, K1 = K1, K2 = K2, e11 = e11, e12 = e12, e21 = e21, e22 = e22, m1 = m1, m2 = m2)

## Attack rate
A <- 2
aii <- A # this will be multiplied by wii above. This helps maintain consistency with previous figure 3.

## confirmed that I don't need to adjusted A (total attack rate investment)
# wii*A + wij*A = (wii + wij)*A = 1*A = A

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

# note though that aij already contains wij implicitly.
# also aii goes up to A/wii because I have to multiply it by wii
concaveDN_df <- data.frame(n = n_concaveDN, aii = seq(0,A/wii,0.01)) %>% mutate(aij = A*(1-(aii*wii/A)^n_concaveDN)^(1/n_concaveDN))
linear_df <- data.frame(n = n_linear, aii = seq(0,A/wii,0.01)) %>% mutate(aij = A*(1-(aii*wii/A)^n_linear)^(1/n_linear))
concaveUP_df <- data.frame(n = n_concaveUP, aii = seq(0,A/wii,0.01)) %>% mutate(aij = A*(1-(aii*wii/A)^n_concaveUP)^(1/n_concaveUP))

bind_rows(concaveDN_df, linear_df, concaveUP_df) %>%
  ggplot(., aes(x=aii*wii, y=aij, color = factor(n))) + geom_line() + ggtitle("Attack rate trade-offs")

bind_rows(concaveDN_df, linear_df, concaveUP_df) %>%
  ggplot(., aes(x=aii*wii/(aii*wii+aij), y=aii*wii+aij, color = factor(n))) + geom_line() + ggtitle("Relationship between specialization and total attack rate")

######################################### MACARTHUR MODEL #######################################################

## Mutant growth rates
# note that aijm already include wii, so I keep this equation the same
mC1_MacArthur <- expression(e11*a11m*R1 + e12*a12m*R2 - m1)
mC2_MacArthur <- expression(e21*a21m*R1 + e22*a22m*R2 - m2)

#### PARAMETERS ####

# adding wii here. Interestingly, there is no opportunity to add
# wij. This suggests that the tradeoff between wii and wij is effectively removed,
# and wij simply adjusts the attack rate. In this case, since aii is already A*0.61
# multiplying it by wii=0.6, makes a11 = a22 = 0.7. I think this should flip
# specialization the other way

# using Mathematica, I identified an alternative tradeoff function that includes wij

# trying wij

## Linear
params_2C_2R_MacArthur_linear <- data.frame(general_parameters,
                                            a11 = aii*wii, 
                                            a12 = eval(a12_tradeoff, data.frame(A=A, a11m=aii*wii, n=n_linear)),                     
                                            a21 = eval(a21_tradeoff, data.frame(A=A, a22m=aii*wii, n=n_linear)), 
                                            a22 = aii*wii)                      
params_1C_2R_MacArthur_linear <- params_2C_2R_MacArthur_linear[c("r1","r2","K1","K2","e11","e12","a11","a12","m1")]


## concaveDN
params_2C_2R_MacArthur_concaveDN <- params_2C_2R_MacArthur_linear
params_2C_2R_MacArthur_concaveDN["a12"] <- eval(a12_tradeoff, data.frame(A=A, a11m=aii*wii, n=n_concaveDN))
params_2C_2R_MacArthur_concaveDN["a21"] <- eval(a21_tradeoff, data.frame(A=A, a22m=aii*wii, n=n_concaveDN))
params_1C_2R_MacArthur_concaveDN <- params_2C_2R_MacArthur_concaveDN[c("r1","r2","K1","K2","e11","e12","a11","a12","m1")]

## concaveUP
params_2C_2R_MacArthur_concaveUP <- params_2C_2R_MacArthur_linear
params_2C_2R_MacArthur_concaveUP["a12"] <- eval(a12_tradeoff, data.frame(A=A, a11m=aii*wii, n=n_concaveUP))
params_2C_2R_MacArthur_concaveUP["a21"] <- eval(a21_tradeoff, data.frame(A=A, a22m=aii*wii, n=n_concaveUP))
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

## Gather data
sim_data <- bind_rows(
  mutate(sim_2C_2R_MacArthur_linear, Model = "MacArthur", n = n_linear, Competitor = "Yes", a_eff_C1 = a11+a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_MacArthur_linear, Model = "MacArthur", n = n_linear, Competitor = "No", a_eff_C1 = a11+a12, Ctot = C1),
  mutate(sim_2C_2R_MacArthur_concaveDN, Model = "MacArthur", n = n_concaveDN, Competitor = "Yes", a_eff_C1 = a11+a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_MacArthur_concaveDN, Model = "MacArthur", n = n_concaveDN, Competitor = "No", a_eff_C1 = a11+a12, Ctot = C1),
  mutate(sim_2C_2R_MacArthur_concaveUP, Model = "MacArthur", n = n_concaveUP, Competitor = "Yes", a_eff_C1 = a11+a12, Ctot = C1+C2), 
  mutate(sim_1C_2R_MacArthur_concaveUP, Model = "MacArthur", n = n_concaveUP, Competitor = "No", a_eff_C1 = a11+a12, Ctot = C1))

#### PREDICTING EFFECTS OF CHARACTER DISPLACEMENT ON C-R DYNAMICS ####

tradeoff_df <- bind_rows(mutate(concaveDN_df, n = n_concaveDN), 
                         mutate(linear_df, n = n_linear), 
                         mutate(concaveUP_df, n = n_concaveUP)) %>%
  data.frame(., general_parameters) %>%
  # note I only multiply aii by wii, because aij implicitly has wij from tradeoff function
  mutate(a11 = aii*wii, a12 = aij, a21 = aij, a22 = aii*wii)

## MacArthur
MacArthur_displacements <- bind_rows(
  mutate(sim_2C_2R_MacArthur_linear, Model = "MacArthur", n = n_linear, Competitor = "Yes")[c(1,nrow(sim_2C_2R_MacArthur_linear)), ], 
  mutate(sim_1C_2R_MacArthur_linear, Model = "MacArthur", n = n_linear, Competitor = "No")[c(1,nrow(sim_1C_2R_MacArthur_linear)), ],
  mutate(sim_2C_2R_MacArthur_concaveDN, Model = "MacArthur", n = n_concaveDN, Competitor = "Yes")[c(1,nrow(sim_2C_2R_MacArthur_concaveDN)), ], 
  mutate(sim_1C_2R_MacArthur_concaveDN, Model = "MacArthur", n = n_concaveDN, Competitor = "No")[c(1,nrow(sim_1C_2R_MacArthur_concaveDN)), ],
  mutate(sim_2C_2R_MacArthur_concaveUP, Model = "MacArthur", n = n_concaveUP, Competitor = "Yes")[c(1,nrow(sim_2C_2R_MacArthur_concaveUP)), ], 
  mutate(sim_1C_2R_MacArthur_concaveUP, Model = "MacArthur", n = n_concaveUP, Competitor = "No")[c(1,nrow(sim_1C_2R_MacArthur_concaveUP)), ]) 
MacArthur_displacements$Time <- rep(c("Begin","End"),6)

MacArthur_2C_2R_tradeoff_dynamics <- CR_dynamics(init_parameters = select(tradeoff_df, r1:a22), init_states = states_2C_2R, eco_CR_model = MacArthur_2C_2R)
MacArthur_1C_2R_tradeoff_dynamics <- CR_dynamics(init_parameters = select(tradeoff_df, r1:a12), init_states = states_1C_2R, eco_CR_model = MacArthur_1C_2R)

## Plot

# Note that I only adjust aii to wii*aii in the tradeoff_df. This ensures that aii and aij are on the same scale
# since aij already incorporates wij implicitly (see previous trade-offs section)
# I don't adjust a11 or 12. This is because I already adjusted a11 in the eco-evo simulation
# and a12 is adjusted automatically (it's determined by a11)

# attack rates
Mac_attack <- ggplot(tradeoff_df, aes(x=wii*aii/(wii*aii+aij), y=wii*aii+aij, color = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=a11+a12, fill=factor(n), shape=Competitor), size=3, color="white") + 
  ylab(expression(a[ii]+a[ij])) +
  scale_shape_manual(values=c(24,21)) +
  scale_color_manual(name="Tradeoff (n)", values=cbPalette)  +
  scale_fill_manual(name="Tradeoff (n)", values=cbPalette) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

# resources
Mac_resources <- cbind.data.frame(select(tradeoff_df, n, aii, aij), MacArthur_2C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=wii*aii/(wii*aii+aij), y=R1+R2, color = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=R1+R2, fill=factor(n), shape=Competitor), size=3, color="white") +
  ylab(expression(R[i]+R[j])) +
  scale_shape_manual(values=c(24,21)) +
  scale_color_manual(name="Tradeoff (n)", values=cbPalette)  +
  scale_fill_manual(name="Tradeoff (n)", values=cbPalette) +
  guides(fill = guide_legend(override.aes=list(shape=21)))


# stability
Stability_Mac <- cbind.data.frame(select(tradeoff_df, n, aii, aij), MacArthur_2C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=wii*aii/(wii*aii+aij), y=-1*max.Re.eigen, color = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=-1*max.Re.eigen, fill=factor(n), shape=Competitor), size=3, color="white", show.legend=FALSE) + 
  geom_point(data = filter(MacArthur_displacements, Time=="Begin"), aes(x = a11/(a11+a12), y=-1*max.Re.eigen, fill=factor(n), shape=Competitor), size=1.5, color="white") +
  ylab(expression("Stability "(-lambda[max]))) + 
  scale_shape_manual(values=c(24,21))  +
  scale_color_manual(name="Tradeoff (n)", values=cbPalette)  +
  scale_fill_manual(name="Tradeoff (n)", values=cbPalette) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

# merge plots
Macs <- plot_grid(
  Mac_attack + xlab("") + ylab(bquote(atop("Effective attack rate","("*a[1*","*j]+a[2*","*j]*")"))) +
    theme_cowplot() +
    theme(legend.position = "none", axis.text = element_text(size=10), axis.text.x = element_blank(), plot.title = element_text(size = 10), axis.title.y=element_text(vjust=0, size=12)), 
  Mac_resources + xlab("") + ylab(bquote(atop("Resource abundance","("*R[1]+R[2]*")"))) + 
    theme_cowplot() +
    theme(legend.position = "none", axis.text = element_text(size=10), axis.text.x = element_blank(), axis.title.y=element_text(vjust=0, size=12)), 
  Stability_Mac + scale_x_continuous(labels = c(0,0.25,0.5,0.75,1)) + xlab("") + ylab(bquote(atop("Food-web stability","("*-lambda[max]*")"))) + 
    theme_cowplot() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(legend.position = "none", axis.text = element_text(size=10), axis.title.y=element_text(vjust=0, size=12)), 
  ncol=1, align='v', labels = c("d)","e)","f)"), hjust = -4, vjust = 0.5, label_size = 10
)
Macs_title <- ggdraw() + draw_label("Resources: Different Habitat", fontface = 'bold', size = 10, hjust = 0.2)
plot_Macs <- plot_grid(Macs_title, Macs, ncol = 1, rel_heights = c(0.1, 1))

plot_Macs
#save_plot(filename="Fig_3_MacArthur_aiiwii_response_to_reviewer.pdf", plot=plot_Macs, base_height = 6.5, base_width = 6)


# check three species situation
# this makes it difficult for the normalization to "work"
cbind.data.frame(select(tradeoff_df, n, aii, aij), MacArthur_1C_2R_tradeoff_dynamics) %>%
  ggplot(., aes(x=wii*aii/(wii*aii+aij), y=R1+R2, color = factor(n))) + 
  geom_line() +
  geom_point(data = filter(MacArthur_displacements, Time=="End"), aes(x = a11/(a11+a12), y=R1+R2, fill=factor(n), shape=Competitor), size=3, color="white") +
  ylab(expression(R[i]+R[j])) +
  scale_shape_manual(values=c(24,21)) +
  scale_color_manual(name="Tradeoff (n)", values=cbPalette)  +
  scale_fill_manual(name="Tradeoff (n)", values=cbPalette) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

# check what determines zero eigenvalues

min_eigen <- min(-1*MacArthur_2C_2R_tradeoff_dynamics$max.Re.eigen)

cbind.data.frame(select(tradeoff_df, n, aii, aij), MacArthur_2C_2R_tradeoff_dynamics) %>%
  filter(-1*max.Re.eigen == min_eigen)

# note that a11 and a12 as well as a21 and a22 are virtually identical.
filter(MacArthur_2C_2R_tradeoff_dynamics, a11 < 1.02 & a11 > 0.98)

