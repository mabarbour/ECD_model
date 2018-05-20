
#### A general function for an eco-evo simulation using Adaptive Dynamics ####

eco_evo_CR <- function(init_parameters,              # data frame of initial parameter values.
                       init_states,                  # data frame of initial state variables
                       eco_CR_model,                 # named function that describes the ecological Consumer-Resource Dynamics
                       AD_CR_models,                 # vector of named functions describing the Adaptive Dynamics of the Consumer-Resource Model
                       species_traits,               # vector that matches each trait to the mutant dynamics in AD_CR_models 
                       eco_param_names,              # vector of evolving parameter names in the ecological model.
                       eco_tradeoff_names,           # vector of trade-off parameter names in the ecological model. Order must correspond to 'eco_param_names'
                       evo_param_names,              # vector of mutant parameter names in the Adaptive Dynamics model. Order must correspond to 'eco_param_names'
                       evo_tradeoff_names,           # vector of mutant trade-off parameter names in the Adaptive Dynamics model. Order must correspon to 'eco_param_names'
                       tradeoff_functions,           # names function of the trade-off with corresponding parameters
                       extra_tradeoff_params=c(NOT_USED=1), # vector of extra tradeoff parameters. Default is a nonsense parameter to just make sure the function runs
                       mut.size=0.001,               # absolute change in evolving parameter. Directionality is chosen at random.
                       mut.suc.thres=1e-4,           # threshold of mutant success rate for identifying ESS.
                       j.step=500,                   # number of evolutionary time steps over which to assess evolutionary change.
                       abund.thres=1e-4,             # abundance threshold for feasible equilibriums
                       mut.reps=1000) {             # number of mutation iterations to try for simulation. Note that this is the upper limit of the simulation length. Usually the traits will hit their ESS or a constraint on the parameter space before this limit is reached.   
                            
  #### SETUP ####
  
  ## Get steady-state resource and consumer abundances from initial parameter values and state variables 
  init.df <- cbind.data.frame(init_parameters, init_states)  # create data frame of initial parameters and states
  set.temp <- as.matrix(init.df[1, ])                        # convert to format for 'safe.runsteady' function
  set.steady <- safe.runsteady(set.temp[1,colnames(init_states)], func = eco_CR_model, parms = set.temp[1,colnames(init_parameters)]) # get steady-state
  eq.jac <- jacobian.full(y = set.steady$y, func = eco_CR_model, parms = set.temp[1,colnames(init_parameters)]) # evaluate local stability at steady-state
  
  ## Update with new steady state and eigenvalues after the initial run
  init.df[ ,colnames(init_states)] <- set.steady$y
  init.df$max.Re.eigen <- max(Re(eigen(eq.jac)$values))
  init.df$max.Im.eigen <- max(Im(eigen(eq.jac)$values))
  
  
  #### ECO-EVO SIMULATION ####
  
  df <- init.df            # start with initial df
  sim.list <- list()       # generate list for simulations to fill in
  for(i in 1:dim(df)[1]){
    
    set <- as.matrix(df[i, ])  # choose parameter set to examine 
    sub.sim.list <- list()     # setup and run sub-simulations
    
    ## Input first data vector
    sub.sim.list[[1]] <- c(set[1, ], 
                           feasibility = 1, 
                           steady.state = 1, 
                           species = NA, 
                           mut.suc = NA, 
                           p.try = NA, 
                           sim.number = i, 
                           sequence = 1) 
    
    ## Run simulation
    for(j in 2:mut.reps){
      
      test <- sub.sim.list %>% plyr::ldply()   # get data frame of output to test whether to keep the simulation running or stop it
      
      ## If the prior simulation results in the exclusion of a species, stop the simulation
      if(sub.sim.list[[j-1]]["feasibility"] == 0) {
        print("Species went extinct!")
        break
      }
      
      ## If the prior simulation results in an unstable system, stop the simulation
      if(sub.sim.list[[j-1]]["max.Re.eigen"] > 0) {
        print("System no longer stable!")
        break
      }
      
      ## If there are no successful mutants after a defined number of mutation steps (j.step), stop the simulation
      if(j < j.step){
        mut.suc.rate <- 1 
      }
      if(j > j.step){
        mut.suc.rate <- mean(test[(j-j.step):(j-1),"mut.suc"], na.rm = TRUE)
      } 
      #if(is.na(mut.suc.rate)==TRUE){
      #  print("what is going on?")
      #}
      if(is.na(mut.suc.rate) & any(test[nrow(test),c(eco_param_names, eco_tradeoff_names)] <= mut.size) == TRUE) { # colnames(init_parameters)
        #sub.sim.list[[j]] <- sub.sim.list[[j-1]]
        #next
        #if(j < mut.reps){
        #  next
        #}
        #if(j == mut.reps){
        #  print("Constraint reached!")
        #}
        print("Constraint reached!") 
        break 
      }
      if(all(mut.suc.rate < mut.suc.thres) & all(test[nrow(test),c(eco_param_names, eco_tradeoff_names)] > mut.size) == TRUE) { # colnames(init_parameters)
        #sub.sim.list[[j]] <- sub.sim.list[[j-1]]
        #if(j < mut.reps){
        #  next
        #}
        #if(j == mut.reps){
        #  print("ESS reached!")
        #}
        #next
        print("ESS reached!") 
        break 
      }
      if(j == mut.reps) {
        print("Maximum mutation iterations reached!")
        break
      }
      
      ## If the prior simulation is at a feasible and steady state, proceed with simulation
      if(sub.sim.list[[j-1]]["feasibility"] == 1 & 
         sub.sim.list[[j-1]]["max.Re.eigen"] < 0){
        
        ## Extract parameters and state variables
        parameters <- sub.sim.list[[j-1]][colnames(init_parameters)]
        state <- sub.sim.list[[j-1]][colnames(init_states)]
        
        ## Setup traits for testing mutant growth rates 
        sp <- sample(1:length(AD_CR_models), 1)                                                 # randomly choose which mutant to evaluate its growth rate
        random.trait <- sample(1:length(eco_param_names), 1)                                    # randomly sample one trait to evolve
        eco.trait.name <- eco_param_names[random.trait]                                         # get the name of the evolving ecological trait
        eco.tradeoff.name <- eco_tradeoff_names[random.trait]                                   # get the name of the trade-off trait
        eco.other.name <- c(eco_param_names[-random.trait], eco_tradeoff_names[-random.trait])  # get the names of the non-evolving traits. Important for calculating mutant growth rates.
        evo.trait.name <- evo_param_names[random.trait]                                         # get the name of the mutant trait
        evo.tradeoff.name <- evo_tradeoff_names[random.trait]                                   # get the name of the mutant trade-off trait
        evo.other.name <- c(evo_param_names[-random.trait], evo_tradeoff_names[-random.trait])  # get the names of the non-evolving mutant traits. Important for calculating mutant growth rates.
        mut <- ifelse(rbinom(1,1,0.5) == 0, -mut.size, mut.size)                                # randomly determine the direction of the mutation (positive or negative change)
        
        ## Get mutant trait values and calculate its growth rate assuming it is rare in the popuplation 
        new.ps <- parameters
        new.ps[evo.trait.name] <- new.ps[eco.trait.name] + mut                                         # get mutant trait value
        new.ps[evo.tradeoff.name] <- eval(expr = tradeoff_functions[random.trait], #[[1]],                       # get mutant trade-off trait value
                                          envir = data.frame(t(c(new.ps, mut=mut, extra_tradeoff_params))))
        new.ps[evo.other.name] <- new.ps[eco.other.name]                                               # keep non-evolving traits the same for the mutant
        eigen.mut <- eval(AD_CR_models[species_traits[random.trait]], data.frame(t(c(new.ps, state)))) # calculate mutant growth rate
        
        ## If the mutant can invade and all traits are positive, update the ecological trait values and find the new equilibrium abundances
        if(eigen.mut > 0 & all(new.ps[c(eco_param_names, eco_tradeoff_names)] >= mut.size) == TRUE){ # trying mut.size instead of zero, because there is some issues with NaN values at very small values
          
          parameters[eco.trait.name] <- new.ps[evo.trait.name] 
          parameters[eco.tradeoff.name] <- new.ps[evo.tradeoff.name]

          ## Run simulation with updated parameters
          out <- safe.runsteady(y = state, func = eco_CR_model, parms = parameters)
          
          ## If all species persist, update the simulation output with the evolved parameters and details on the stability of the system
          if(any(out$y < abund.thres) == FALSE){
            
            eq.jac <- jacobian.full(y = out$y, func = eco_CR_model, parms = parameters) 
            
            sub.sim.list[[j]] <- c(parameters, out$y,  
                                     max.Re.eigen = max(Re(eigen(eq.jac)$values)), 
                                     max.Im.eigen = max(Im(eigen(eq.jac)$values)),
                                     feasibility = 1,
                                     steady.state = attr(out, "steady"),
                                     species = sp, mut.suc = 1, p.try = 1,
                                     sim.number = i, sequence = j)
          }
          
          ## If any species goes extinct, set eigenvalues and steady-state to NA because I'm not interested in these values when a species has been excluded.
          if(any(out$y < abund.thres) == TRUE){
            
            sub.sim.list[[j]] <- c(parameters, out$y, 
                                     max.Re.eigen = NA,
                                     max.Im.eigen = NA,
                                     feasibility = 0, 
                                     steady.state = NA,
                                     species = sp, mut.suc = 1, p.try = 1,
                                     sim.number = i, sequence = j)
          }
        }
        
        ## If the mutant can't invade and all trait values are positive, keep the trait, equilibrium, and stability values the same as the previous iteration
        if(eigen.mut < 0 & all(new.ps[c(eco_param_names, eco_tradeoff_names)] >= mut.size) == TRUE){ 
          
          sub.sim.list[[j]] <- c(sub.sim.list[[j-1]][c(colnames(init_parameters), colnames(init_states), "max.Re.eigen", "max.Im.eigen","feasibility","steady.state")],  
                                   species = sp, mut.suc = 0, p.try = 1, 
                                   sim.number = i, sequence = j)
        }
        
        ## If any of the trait values becomes negative, keep the trait, equilibrium, and stability values the same as the previous iteration
        if(all(new.ps[c(eco_param_names, eco_tradeoff_names)] >= mut.size) == FALSE){ 
          
          sub.sim.list[[j]] <- c(sub.sim.list[[j-1]][c(colnames(init_parameters), colnames(init_states), "max.Re.eigen", "max.Im.eigen","feasibility","steady.state")],  
                                 species = sp, mut.suc = NA, p.try = -1, 
                                 sim.number = i, sequence = j)
        }
      }
    }
    
    ## Turn the sub-simulations into a data frame and output to the simulation list.
    sim.list[[i]] <- sub.sim.list %>% plyr::ldply()
  }
  
  ## Convert the simulation list into a data frame  
  eco.evo.sim.df <- sim.list %>% plyr::ldply() 
  return(eco.evo.sim.df)
  
}