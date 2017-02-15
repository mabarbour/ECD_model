## load required libraries
library(tidyverse)
library(cowplot)

## load required data ----
C1.ESS.df <- read.csv('data/evol.symmetry.df.C1.3sp.csv')

## source required functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

#### Setup the initial state variables and global parameters for all of the models ----

## set state variables for all models
R1 <- 2
R2 <- 2
C1 <- 1
C2 <- 1
i.state.2C_2R <- c(R1 = R1, R2 = R2, C1 = C1, C2 = C2) 

## set duration of simulations
Time <- 1000

## set up general code for bifurcation plot. Maintaining 'w' as the key parameter to manipulate
aii <- C1.ESS.df$a11[1]
aij <- C1.ESS.df$a12[1]

aii_seq <- seq(aii, aii*2, by = 0.01)[1:75]
aij_seq <- seq(aij, 0, by = -0.01)[1:75]

w11 <- w22 <- C1.ESS.df$w11[1]
#a_alt <- a_start
aECD <- w11*aii_seq + (1-w11)*aij_seq

## Bifurcation plot code
# this code loops over the sequence of different parameter values
out <- list()
for (i in 1:length(aii_seq)) {
  # general parameters. Taken from Figure 3 in McCann et al. 2005. 
  r1 <- r2 <- 1
  K1 <- K2 <- 3.25
  e11 <- e12 <- e21 <- e22 <- 0.8
  h11 <- h12 <- h21 <- h22 <- 0.4
  m1 <- m2 <- 1
  
  # parameters related to aECD
  #w <- 0.8
  #a_prime <- param.seq.a_prime[i]
  #a_alt <- 2
  
  # set params 
  param.loop <- c(r1 = r1, r2 = r2, K1 = K1, K2 = K2,
                  e11 = e11, e12 = e12, e21 = e21, e22 = e22,
                  h11 = h11, h12 = h12, h21 = h21, h22 = h22,
                  a11 = aii_seq[i], a12 = aij_seq[i], a21 = aij_seq[i], a22 = aii_seq[i],
                  m1 = m1, m2 = m2, w11 = w11, w22 = w22)
  
  # Run the experiment.
  init <- ode(i.state.2C_2R, 
              1:Time, ECD_model, param.loop)
  # Rerun the experiment with the state variables at the end of the initial simulation. This will enable me to determine the range of state variables around the equilibrium.
  out[[i]] <- ode(init[Time,-1], 1:Time, ECD_model, param.loop)[,-1]   
}

sim.list <- list()
for (i in 1:length(aii_seq)) {
  sim.list[[i]] <- data.frame(rep = i,
                              C1_max = max(out[[i]][900:1000,"C1"]),
                              C1_min = min(out[[i]][900:1000,"C1"]),
                              R1_max = max(out[[i]][900:1000,"R1"]),
                              R1_min = min(out[[i]][900:1000,"R1"]))
}


sim.df <- ldply(sim.list) %>%
  mutate(aECD = aECD) %>%
  gather(rep, aECD)
colnames(sim.df) <- c("rep","aECD","temp","value")
sim.df <- separate(sim.df, temp, into = c("species","max_min"))

## Plot results from the simulation
bifur <- ggplot(sim.df, aes(x = rep, y = value, color = species)) +
  geom_point(size = 1) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), guide = "none") +
  labs(y = "Population density (max and min)",
       x = "Evolutionary time steps") +
  theme_cowplot(font_size = 18)

#save_plot("amnat_bifur_plot.pdf", bifur, base_height = 8.5, base_width = 11)

