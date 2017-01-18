#### Description: 4 species model with only one consumer (C1) evolving simulation

## source general parameters and functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

## load data ----
all.feas.stable.df <- read.csv('data/all.feas.stable.df.csv') %>% 
  tbl_df() %>%
  select(r1:C2, max.Re.eigen, max.Im.eigen)


## setup and run simulation ----
df <- all.feas.stable.df
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
    if(j < 100){
      delta.a11 <- 1 # arbitrary number greater than threshold value
      delta.a22 <- 1 # arbitrary number greater than threshold value
    }
    if(j > 100){
      delta.a11 <- mean(diff(test[(j-100):(j-1),"a11"]))
      delta.a22 <- mean(diff(test[(j-100):(j-1),"a22"]))
    } 
    if(delta.a11 < trait.thres & delta.a22 < trait.thres) break 
    
    # if the prior simulation is at a feasible and steady state, proceed with simulation
    if(sub.sim.list.4[[j-1]]["feas.4sp"] == 1 & 
       sub.sim.list.4[[j-1]]["max.Re.eigen"] < 0){
      
      # extract parameters and state variables
      parameters <- sub.sim.list.4[[j-1]][1:20]
      state <- sub.sim.list.4[[j-1]][c("R1","R2","C1","C2")]
      
      # randomly assign mutation in attack rate to one of the species
      sp <- 1 # only allow C1 to evolve
      mut <- ifelse(rbinom(1,1,0.5) == 0, -mut.size, mut.size)
      
      # determine whether mutant can invade or not by calculating whether it has positive growth rate (i.e. positive eigen value)
      new.ps <- parameters
      
      if(sp == 1){
        new.ps["a11m"] <- new.ps["a11"] + mut
        new.ps["a12m"] <- new.ps["a12"] - mut
        
        jac.mut <- jacobian.full(y = c(state, mC1 = 0), func = mutC1_ECD_model, parms = new.ps)
        eigen.mut <- max(Re(eigen(jac.mut)$values)) 
      }
      
      # if the mutant can invade and attack rates are positive, update attack rates and find the new equilibrium abundances
      if(eigen.mut > 0 & all(new.ps >= 0) == TRUE){
        if(sp == 1){
          parameters["a11"] <- new.ps["a11m"]
          parameters["a12"] <- new.ps["a12m"]
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


## convert list into data frame then write the data to a new file ----
evol.sim.df.4.1evo <- sim.list.4 %>% ldply() 

write.csv(x = evol.sim.df.4.1evo, file = 'data/evol.sim.df.4sp.1evo.csv')

