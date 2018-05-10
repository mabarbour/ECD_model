CR_dynamics <- function(init_parameters,              # data frame of initial parameter values.
                        init_states,                  # data frame of initial state variables
                        eco_CR_model,                 # named function that describes the ecological Consumer-Resource Dynamics
                        abund.thres=1e-4) { # 1e-2             # abundance threshold for feasible equilibriums
  
  init.df <- data.frame(init_parameters, init_states, max.Re.eigen = NA, max.Im.eigen = NA)  # create data frame of initial parameters and states
  
  #out.list <- list()
  for(i in 1:nrow(init.df)){
    
    ## Get steady-state resource and consumer abundances from initial parameter values and state variables 
    set.temp <- as.matrix(init.df[i, ])                        # convert to format for 'safe.runsteady' function
    set.steady <- safe.runsteady(set.temp[1,colnames(init_states)], func = eco_CR_model, parms = set.temp[1,colnames(init_parameters)]) # get steady-state
    eq.jac <- jacobian.full(y = set.steady$y, func = eco_CR_model, parms = set.temp[1,colnames(init_parameters)]) # evaluate local stability at steady-state
    
    ## Update with new steady state and eigenvalues after the initial run
    init.df[i,colnames(init_states)] <- set.steady$y
    init.df[i,"max.Re.eigen"] <- max(Re(eigen(eq.jac)$values))
    init.df[i,"max.Im.eigen"] <- max(Im(eigen(eq.jac)$values))
  }
  CR_dynamics_output <- init.df
  return(CR_dynamics_output)
}