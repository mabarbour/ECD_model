
#### Description: 4 species model with consumer coevolution simulation,
####              but this is only for one parameter set where resource
####              parameters are equivalent and consumers are mirror
####              images of each other.

## source general parameters and functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

trait.thres <- 1e-10
j.step <- 1000 # 100 for other simulations

## general parameters for simulations
r <- 1
K <- 4
e <- 0.8
h <- 0.4
aii <- 2
aij <- 1.2
#a <- 2
w <- 0.6
m <- 1
R <- 2
C <- 1


## 4 species simulation ----

## create initial data
sym4.init.df <- data.frame(r1 = r, r2 = r, K1 = K, K2 = K,
                         e11 = e, e12 = e, e21 = e, e22 = e,
                         h11 = h, h12 = h, h21 = h, h22 = h,
                         a11 = aii, a12 = aij, a21 = aij, a22 = aii,
                         w11 = w, w22 = w, m1 = m, m2 = m,
                         R1 = R, R2 = R, C1 = C, C2 = C)

set.temp <- as.matrix(sym4.init.df[1, ])

set.steady <- safe.runsteady(set.temp[1,c("R1","R2","C1","C2")], func = ECD_model, parms = set.temp[1,1:20])

sym4.eq.jac <- jacobian.full(y = set.steady$y, 
                           func = ECD_model, parms = set.temp[1,1:20])

## update with new steady state and eigenvalues after the initial run
sym4.init.df[ ,c("R1","R2","C1","C2")] <- set.steady$y
sym4.init.df$max.Re.eigen <- max(Re(eigen(sym4.eq.jac)$values))
sym4.init.df$max.Im.eigen <- max(Im(eigen(sym4.eq.jac)$values))


## setup and run simulation
df <- sym4.init.df
sim.list.4 <- list()
for(i in 1:dim(df)[1]){
  
  # choose parameter set to examine and calculate the maximum real eigenvalue
  set <- as.matrix(df[i, ])
  
  ## setup and run sub-simulations
  sub.sim.list.4 <- list()
  
  # input first data vector
  sub.sim.list.4[[1]] <- c(set[1, ], feas.4sp = 1, steady.4sp = 1, 
                           species = NA, mut.suc = NA, p.try = NA, 
                           sim.number = i, sequence = 1) 
  
  # simulation
  for(j in 2:mut.reps){
    
    test <- sub.sim.list.4 %>% ldply()
    
    # if the prior simulation results in the exclusion of a species, stop the simulation
    if(sub.sim.list.4[[j-1]]["feas.4sp"] == 0) break 
    
    # if the prior simulation results in an unstable system, stop the simulation
    if(sub.sim.list.4[[j-1]]["max.Re.eigen"] > 0) break  
    
    # if we are no longer observing substantial trait change after 100 mutations step, stop the simulation
    if(j < j.step){
      delta.a11 <- 1 # arbitrary number greater than threshold value
      delta.a22 <- 1 # arbitrary number greater than threshold value
    }
    if(j > j.step){
      delta.a11 <- mean(diff(test[(j-j.step):(j-1),"a11"]))
      delta.a22 <- mean(diff(test[(j-j.step):(j-1),"a22"]))
    } 
    if(delta.a11 < trait.thres & delta.a22 < trait.thres) break 
    
    # if the prior simulation is at a feasible and steady state, proceed with simulation
    if(sub.sim.list.4[[j-1]]["feas.4sp"] == 1 & 
       sub.sim.list.4[[j-1]]["max.Re.eigen"] < 0){
      
      # extract parameters and state variables
      parameters <- sub.sim.list.4[[j-1]][1:20]
      state <- sub.sim.list.4[[j-1]][c("R1","R2","C1","C2")]
      
      # randomly assign mutation in attack rate to one of the species
      sp <- rbinom(1,1,0.5) + 1
      mut <- ifelse(rbinom(1,1,0.5) == 0, -mut.size, mut.size)
      
      # determine whether mutant can invade or not by calculating whether it has positive growth rate (i.e. positive eigen value)
      new.ps <- parameters
      
      if(sp == 1){
        new.ps["a11m"] <- new.ps["a11"] + mut
        new.ps["a12m"] <- new.ps["a12"] - mut
        
        jac.mut <- jacobian.full(y = c(state, mC1 = 0), func = mutC1_ECD_model, parms = new.ps)
        eigen.mut <- max(Re(eigen(jac.mut)$values)) 
      }
      if(sp == 2){
        new.ps["a22m"] <- new.ps["a22"] + mut
        new.ps["a21m"] <- new.ps["a21"] - mut
        
        jac.mut <- jacobian.full(y = c(state, mC2 = 0), func = mutC2_ECD_model, parms = new.ps)
        eigen.mut <- max(Re(eigen(jac.mut)$values)) # jac.mut[6,6] # 
      }
      
      # if the mutant can invade and attack rates are positive, update attack rates and find the new equilibrium abundances
      if(eigen.mut > 0 & all(new.ps >= 0) == TRUE){
        if(sp == 1){
          parameters["a11"] <- new.ps["a11m"]
          parameters["a12"] <- new.ps["a12m"]
        }
        if(sp == 2){
          parameters["a21"] <- new.ps["a21m"]
          parameters["a22"] <- new.ps["a22m"]
        }
        
        # run simulation with new parameters
        out <- safe.runsteady(y = state, func = ECD_model, parms = parameters)
        
        # return steady states and maximum real eigenvalues
        if(any(out$y < abund.thres) == FALSE){
          
          eq.jac <- jacobian.full(y = out$y, func = ECD_model, parms = parameters) 
          
          sub.sim.list.4[[j]] <- c(parameters, out$y,  
                                   max.Re.eigen = max(Re(eigen(eq.jac)$values)), 
                                   max.Im.eigen = max(Im(eigen(eq.jac)$values)),
                                   feas.4sp = 1,
                                   steady.4sp = attr(out, "steady"),
                                   species = sp, mut.suc = 1, p.try = 1,
                                   sim.number = i, sequence = j)
        }
        if(any(out$y < abund.thres) == TRUE){
          
          # set eigenvalues and steady-state to NA because I'm not interested in these values when a species has been excluded.
          sub.sim.list.4[[j]] <- c(parameters, out$y, 
                                   max.Re.eigen = NA,
                                   max.Im.eigen = NA,
                                   feas.4sp = 0, 
                                   steady.4sp = NA,
                                   species = sp, mut.suc = 1, p.try = 1,
                                   sim.number = i, sequence = j)
        }
      }
      
      # if the mutant can't invade or the evolved attack rates are negative, keep attack rates the same and print the same equilibrium values with maximum real eigenvalue
      if(eigen.mut < 0 & all(new.ps >= 0) == TRUE){
        
        sub.sim.list.4[[j]] <- c(sub.sim.list.4[[j-1]][1:28], 
                                 species = sp, mut.suc = 0, p.try = 1, 
                                 sim.number = i, sequence = j)
      }
      if(all(new.ps >= 0) == FALSE){
        
        sub.sim.list.4[[j]] <- c(sub.sim.list.4[[j-1]][1:28], 
                                 species = sp, mut.suc = NA, p.try = -1, 
                                 sim.number = i, sequence = j)
      }
    }
  }
  # turn sub-simulation into data frame and output to simulation list.
  sim.list.4[[i]] <- sub.sim.list.4 %>% ldply()
}


## convert list into data frame then write the data to a new file 
evol.sim.df.4 <- sim.list.4 %>% ldply() 

write.csv(x = evol.sim.df.4, file = 'data/evol.symmetry.df.4sp.csv')



## 3 species simulation ----

## create initial data
C1.init.df <- data.frame(r1 = r, r2 = r, K1 = K, K2 = K,
                         e11 = e, e12 = e, h11 = h, h12 = h,
                         a11 = aii, a12 = aij, m1 = m, w11 = w, 
                         R1 = R, R2 = R, C1 = C)

set.temp <- as.matrix(C1.init.df[1, ])

set.steady <- safe.runsteady(set.temp[1,c("R1","R2","C1")], func = ECD_model.C1.3sp, parms = set.temp[1,1:12])

C1.eq.jac <- jacobian.full(y = set.steady$y, 
                           func = ECD_model.C1.3sp, parms = set.temp[1,1:12])

C1.init.df[ ,c("R1","R2","C1")] <- set.steady$y
C1.init.df$C1.max.eigen <- max(Re(eigen(C1.eq.jac)$values))

## setup and run simulation
df <- C1.init.df
sim.list.C1 <- list()
for(i in 1:dim(df)[1]){
  
  ## choose parameter set to examine and calculate the maximum real eigenvalue
  set <- as.matrix(df[i, ])
  
  ## setup and run sub-simulations
  sub.sim.list.3 <- list()
  
  # input first data vector
  sub.sim.list.3[[1]] <- c(set[1, ], feas.3sp = 1, steady.3sp = 1, 
                           species = NA, mut.suc = NA, p.try = NA,
                           sim.number = i, sequence = 1) 
  
  # run simulation
  for(j in 2:mut.reps){
    
    test <- sub.sim.list.3 %>% ldply()
    
    # if the prior simulation results in the exclusion of a species, stop the simulation
    if(sub.sim.list.3[[j-1]]["feas.3sp"] == 0) break 
    
    # if the prior simulation results in an unstable system, stop the simulation
    if(sub.sim.list.3[[j-1]]["C1.max.eigen"] > 0) break  
    
    # if we are no longer observing substantial trait change after 100 mutations step, stop the simulation
    if(j < j.step){
      delta.a11 <- 1 # arbitrary number greater than threshold value
    }
    if(j > j.step){
      delta.a11 <- mean(diff(test[(j-j.step):(j-1),"a11"]))
    } 
    if(delta.a11 < trait.thres) break 
    
    # if the prior simulation is at a feasible and steady state, proceed with simulation
    if(sub.sim.list.3[[j-1]]["feas.3sp"] == 1 & 
       sub.sim.list.3[[j-1]]["C1.max.eigen"] < 0){
      
      # extract parameters and state variables
      parameters <- sub.sim.list.3[[j-1]][1:12]
      state <- sub.sim.list.3[[j-1]][c("R1","R2","C1")]
      
      # assign mutation in attack rate to one of the species
      sp <- 1 # assign to C1
      mut <- ifelse(rbinom(1,1,0.5) == 0, -mut.size, mut.size)
      
      # determine whether mutant can invade or not by calculating whether it has positive growth rate (i.e. positive eigen value)
      new.ps <- parameters
      
      if(sp == 1){
        new.ps["a11m"] <- new.ps["a11"] + mut
        new.ps["a12m"] <- new.ps["a12"] - mut
        
        jac.mut <- jacobian.full(y = c(state, mC1 = 0), func = mut_ECD_model.C1.3sp, parms = new.ps) 
        eigen.mut <- max(Re(eigen(jac.mut)$values))
      }
      
      # if the mutant can invade and attack rates are positive, update attack rates and print new equilibrium values
      if(eigen.mut > 0 & all(new.ps >= 0) == TRUE){
        if(sp == 1){
          parameters["a11"] <- new.ps["a11m"]
          parameters["a12"] <- new.ps["a12m"]
        }
        
        # run simulation with new parameters
        out <- safe.runsteady(y = state, func = ECD_model.C1.3sp, parms = parameters)
        
        # return steady states and maximum real eigenvalues
        if(any(out$y < abund.thres) == FALSE){
          
          eq.jac <- jacobian.full(y = out$y, func = ECD_model.C1.3sp, parms = parameters)
          
          sub.sim.list.3[[j]] <- c(parameters, out$y,  
                                   C1.max.eigen = max(Re(eigen(eq.jac)$values)), 
                                   feas.3sp = 1,
                                   steady.3sp = attr(out, "steady"),
                                   species = sp, mut.suc = 1, p.try = 1,
                                   sim.number = i, sequence = j)
        }
        if(any(out$y < abund.thres) == TRUE){
          
          # set eigenvalues and steady-state to NA because I'm not interested in these values when a species has been excluded.
          sub.sim.list.3[[j]] <- c(parameters, out$y, 
                                   C1.max.eigen = NA,
                                   feas.3sp = 0, 
                                   steady.3sp = NA,
                                   species = sp, mut.suc = 1, p.try = 1,
                                   sim.number = i, sequence = j)
        }
      }
      
      # if the mutant can't invade or the evolved attack rates are negative, keep attack rates the same and print the same equilibrium values with maximum real eigenvalue
      if(eigen.mut < 0 & all(new.ps >= 0) == TRUE){
        
        sub.sim.list.3[[j]] <- c(sub.sim.list.3[[j-1]][1:18], 
                                 species = sp, mut.suc = 0, p.try = 1,
                                 sim.number = i, sequence = j)
      }
      if(all(new.ps >= 0) == FALSE){
        
        sub.sim.list.3[[j]] <- c(sub.sim.list.3[[j-1]][1:18], 
                                 species = sp, mut.suc = NA, p.try = -1,
                                 sim.number = i, sequence = j)
      }
    }
  }
  # turn sub-simulation into data frame and output to simulation list.
  sim.list.C1[[i]] <- sub.sim.list.3 %>% ldply()
}

## convert list into data frame, then write data to a new file 
evol.symmetry.df.C1 <- sim.list.C1 %>% ldply()

write.csv(x = evol.symmetry.df.C1, file = 'data/evol.symmetry.df.C1.3sp.csv')




## Plots ----

# Color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# subset data frame to only successful evolutionary time steps
evo.3 <- evol.symmetry.df.C1 %>% #filter(evol.symmetry.df.C1, mut.suc != 0) %>%
  select(sequence, a11, a12, w11, 
         R1, R2, C1,  
         max.Re.eigen = C1.max.eigen) %>%
  mutate(sim.type = "3_sp",
         special.C1R1 = a11*w11/(a11*w11+a12*(1-w11)),
         special.C1R2 = a12*(1-w11)/(a11*w11+a12*(1-w11)))

evo.4 <- evol.sim.df.4 %>% #filter(evol.sim.df.4, mut.suc != 0)
  select(sequence, a11, a12, a21, a22, w11, w22,
         R1, R2, C1, C2, 
         max.Re.eigen) %>%
  mutate(sim.type = "4_sp",
         special.C1R1 = a11*w11/(a11*w11+a12*(1-w11)),
         special.C1R2 = a12*(1-w11)/(a11*w11+a12*(1-w11)),
         special.C2R1 = a21*(1-w22)/(a21*(1-w22)+a22*w22),
         special.C2R2 = a22*w22/(a21*(1-w22)+a22*w22))

sim.dur <- dim(evo.4)[1]
evo.sym.df <- bind_rows(evo.4, evo.3) %>% mutate(sim.type = as.factor(sim.type))

## How does consumer specialization evolve over time?
C1.x <- evo.4$special.C1R1[1] 
C1.y <- evo.4$special.C1R2[1] 
C2.x <- C1.y
C2.y <- C1.x
C1.lab <- paste("italic(C[1])")
C2.lab <- paste("italic(C[2])")

spec <- ggplot(evo.sym.df %>% filter(sequence %in% c(1,sim.dur))) +
  geom_segment(aes(x = special.C1R1[1], xend = special.C1R1[2], y = special.C1R2[1], yend = special.C1R2[2]), 
               arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "open"), color = cbbPalette[6]) + 
  geom_segment(aes(x = special.C2R1[1], xend = special.C2R1[2], y = special.C2R2[1], yend = special.C2R2[2]), 
               arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "open"), color = cbbPalette[3]) +
  geom_segment(aes(x = special.C1R1[3], xend = special.C1R1[4], y = special.C1R2[3], yend = special.C1R2[4]), 
               arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "open"), linetype = "dotted", color = cbbPalette[6]) + 
  geom_point(aes(x = special.C1R1[1], y = special.C1R2[1]), shape = 21, size = 10, fill = "white", color = cbbPalette[6]) +
  geom_point(aes(x = special.C2R1[1], y = special.C2R2[1]), shape = 21, size = 10, fill = "white", color = cbbPalette[3]) +
  annotate("text", label = C1.lab, x = C1.x, y = C1.y, parse = TRUE) +
  annotate("text", label = C2.lab, x = C2.x, y = C2.y, parse = TRUE) +
  #geom_text(aes(x = special.C2R1[1], y = special.C2R2[1]), shape = 21, size = 3, fill = "white") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab(expression(Adaptation~to~italic(R[1]))) +
  ylab(expression(Adaptation~to~italic(R[2])))

## How does the effective attack rate of a consumer evolve over time?
eff <- ggplot(evo.sym.df, aes(x = sequence, group = sim.type, linetype = sim.type)) +
  geom_line(aes(y = a11*w11 + a12*(1-w11)), color = cbbPalette[6], show.legend = FALSE) + # C1
  geom_line(aes(y = a22*w22 + a21*(1-w22)), color = cbbPalette[3], show.legend = FALSE) + # C2
  scale_x_continuous(limits = c(0,sim.dur)) +
  scale_linetype_manual(values = c("dashed","solid")) +
  ylab(expression(Effective~attack~rate~(italic(a[ii]*w[ii]~+~a[ij]*(1-w[ii]))))) +
  xlab("Evolutionary time") # number of mutation attempts.

## How does consumer and resource densities change over evolutionary time? 
tidy.evo.df <- evo.sym.df %>% gather(key = species, value = density, R1:C2)

# only resource densities
res.dens <- ggplot(tidy.evo.df %>% filter(species %in% c("R1","R2")), aes(x = sequence, linetype = sim.type, color = species)) +
  geom_line(aes(y = density), show.legend = FALSE) +
  scale_x_continuous(limits = c(0,sim.dur)) +
  scale_y_continuous(limits = c(0,2.5)) +
  scale_linetype_manual(values = c("dashed","solid")) +
  scale_color_manual(values = c(cbbPalette[6], cbbPalette[3], cbbPalette[4], cbbPalette[2])) +
  ylab("Density at equilibrium") +
  xlab("Evolutionary time") # number of mutation attempts.

# both consumers and resources
ggplot(tidy.evo.df, aes(x = sequence, linetype = sim.type, color = species)) +
  geom_line(aes(y = density)) +
  scale_x_continuous(limits = c(0,sim.dur)) +
  scale_y_continuous(limits = c(0,2.5)) +
  scale_linetype_manual(values = c("dashed","solid")) +
  scale_color_manual(values = c(cbbPalette[6], cbbPalette[3], cbbPalette[4], cbbPalette[2])) +
  ylab("Density at equilibrium") +
  xlab("Evolutionary time") # number of mutation attempts.

# total consumers and total resources
tot.res.dens <- ggplot(evo.sym.df %>% 
         mutate(totC = rowSums(cbind(C1, C2), na.rm = TRUE), 
                totR = R1 + R2) %>% 
         gather(key = species, value = total_density, totC, totR), 
       aes(x = sequence, color = species, linetype = sim.type)) +
  geom_line(aes(y = total_density)) +
  scale_x_continuous(limits = c(0,sim.dur)) +
  scale_y_continuous(limits = c(0,4)) +
  scale_linetype_manual(values = c("dashed","solid")) +
  scale_color_manual(values = c(cbbPalette[6], cbbPalette[4])) +
  ylab("Density at equilibrium") +
  xlab("Evolutionary time") # number of mutation attempts.

# total resources only
tot.res.dens <- ggplot(evo.sym.df %>% 
                         mutate(totR = R1 + R2), 
                       aes(x = sequence, linetype = sim.type)) +
  geom_line(aes(y = totR), show.legend = FALSE) +
  scale_x_continuous(limits = c(0,sim.dur)) +
  scale_y_continuous(limits = c(0,3.5)) +
  scale_linetype_manual(values = c("dashed","solid")) +
  scale_color_manual(values = c(cbbPalette[6], cbbPalette[4])) +
  ylab("Total resource density (at equilibrium)") +
  xlab("Evolutionary time") # number of mutation attempts.

## How does the stability of the system evolve over time?
stab <- ggplot(evo.sym.df, aes(x = sequence, group = sim.type, linetype = sim.type)) +
  geom_line(aes(y = -1*max.Re.eigen), show.legend = FALSE) + 
  scale_x_continuous(limits = c(0,sim.dur)) +
  scale_linetype_manual(values = c("dashed","solid")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab(expression(Stability~(-1*lambda))) +
  xlab("Evolutionary time")

## Integrate into figure
fig_theory <- plot_grid(spec, eff, tot.res.dens, stab, align = "hv", labels = "AUTO")

save_plot("figures/fig_theory.png", fig_theory, base_height = 8.5, base_width = 8.5, base_aspect_ratio = 1)
