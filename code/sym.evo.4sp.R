
#### Description: 4 species model with consumer coevolution simulation,
####              but this is only for one parameter set where resource
####              parameters are equivalent and consumers are mirror
####              images of each other.

## source general parameters and functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

## create initial data ----
sym4.init.df <- data.frame(r1 = 1, r2 = 1, K1 = 3.25, K2 = 3.25,
                         e11 = 0.8, e12 = 0.8, e21 = 0.8, e22 = 0.8,
                         h11 = 0.4, h12 = 0.4, h21 = 0.4, h22 = 0.4,
                         a11 = 2, a12 = 2, a21 = 2, a22 = 2,
                         w11 = 0.8, w22 = 0.8, m1 = 1, m2 = 1,
                         R1 = 2, R2 = 2, C1 = 1, C2 = 1)

set.temp <- as.matrix(sym4.init.df[1, ])

set.steady <- safe.runsteady(set.temp[1,c("R1","R2","C1","C2")], func = ECD_model, parms = set.temp[1,1:20])

sym4.eq.jac <- jacobian.full(y = set.steady$y, 
                           func = ECD_model, parms = set.temp[1,1:20])

## update with new steady state and eigenvalues after the initial run
sym4.init.df[ ,c("R1","R2","C1","C2")] <- set.steady$y
sym4.init.df$max.Re.eigen <- max(Re(eigen(sym4.eq.jac)$values))
sym4.init.df$max.Im.eigen <- max(Im(eigen(sym4.eq.jac)$values))


## setup and run simulation ----
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


## convert list into data frame then write the data to a new file ----
evol.sim.df.4 <- sim.list.4 %>% ldply() 

write.csv(x = evol.sim.df.4, file = 'data/evol.symmetry.df.4sp.csv')