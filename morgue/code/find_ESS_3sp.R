
#### Description: Find ESS for C1 3 sp model when resources have equivalent parameters.

## source required functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

## create initial data ----
C1.init.df <- data.frame(r1 = 1, r2 = 1, K1 = 3.25, K2 = 3.25,
                         e11 = 0.8, e12 = 0.8, h11 = 0.4, h12 = 0.4,
                         a11 = 2, a12 = 2, m1 = 1, w11 = 0.8, 
                         R1 = 2, R2 = 2, C1 = 1)

set.temp <- as.matrix(C1.init.df[1, ])

set.steady <- safe.runsteady(set.temp[1,c("R1","R2","C1")], func = ECD_model.C1.3sp, parms = set.temp[1,1:12])

C1.eq.jac <- jacobian.full(y = set.steady$y, 
                           func = ECD_model.C1.3sp, parms = set.temp[1,1:12])

C1.init.df[ ,c("R1","R2","C1")] <- set.steady$y
C1.init.df$C1.max.eigen <- max(Re(eigen(C1.eq.jac)$values))

## setup and run simulation ----
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
    if(j < 100){
      delta.a11 <- 1 # arbitrary number greater than threshold value
    }
    if(j > 100){
      delta.a11 <- mean(diff(test[(j-100):(j-1),"a11"]))
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

## convert list into data frame, then write data to a new file ----
evol.symmetry.df.C1 <- sim.list.C1 %>% ldply()

write.csv(x = evol.symmetry.df.C1, file = 'data/evol.symmetry.df.C1.3sp.csv')
