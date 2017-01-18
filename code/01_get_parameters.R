
#### Description: Use latin hypercube to generate possible parameters that will result in coexistence equilibriums for 4 species.

## load required functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

## Setup simulation ----

## Specificy the number of parameter sets to generate
n.p <- 5000 # 95 hrs to finish on my computer (almost 4 days)

## Create a function that gets the model parameters and time sequence and outputs the final size of populations
oneRun <- function(R1_0, R2_0, C1_0, C2_0,
                   r1, r2, K1, K2,
                   e11, e12, e21, e22,
                   h11, h12, h21, h22,
                   a11, a12, a21, a22,
                   m1, m2, w11, w22){  
  # give parameters, initial state variables, and simulation length
  parameters <- c(r1 = r1, r2 = r2, K1 = K1, K2 = K2,
                  e11 = e11, e12 = e12, e21 = e21, e22 = e22,
                  h11 = h11, h12 = h12, h21 = h21, h22 = h22,
                  a11 = a11, a12 = a12, a21 = a21, a22 = a22,
                  m1 = m1, m2 = m2, w11 = w11, w22 = w22)
  state <- c(R1 = R1_0, R2 = R2_0, C1 = C1_0, C2 = C2_0)
  times <- c(0, 1000) # length of simulation

  # Use ode integration and index to get final state values. Will check later whether 
  
  
  # using failwith() to make function give output even with a fatal error
  # runsteady() identifies the steady-state of the model.
  #safe.runsteady <- failwith(default = c(R1 = NA, R2 = NA, C1 = NA, C2 = NA, steady = 0), runsteady)
  
  out <- safe.runsteady(y = state, func = ECD_model, parms = parameters)
  #out <- safe.stode(y = state, func = ECD_model, parms = parameters, positive = TRUE) #times = times, 
  
  # return steady states
  if(class(out) == "list"){
    return(c(out$y, steady = attr(out, "steady")))
  }
  if(class(out) == "numeric"){
    return(out)
  }
}

## Create a function with mapply that gets a matrix of parameter combinations (combinations are lines and parameters are columns) and returns the output of the function above for each parameter combination.
modelRun <- function (my.pars) {
  return(mapply(oneRun, 
                my.pars[,1], my.pars[,2], my.pars[,3], my.pars[,4], 
                my.pars[,5], my.pars[,6], my.pars[,7], my.pars[,8],
                my.pars[,9], my.pars[,10], my.pars[,11], my.pars[,12],
                my.pars[,13], my.pars[,14], my.pars[,15], my.pars[,16],
                my.pars[,17], my.pars[,18], my.pars[,19], my.pars[,20],
                my.pars[,21], my.pars[,22], my.pars[,23], my.pars[,24]))
}

## Now prepare the hypercube.

# You have to specify a character vector to name the state variables and the model parameters
factors <- c("R1_0","R2_0","C1_0","C2_0",
             "r1","r2","K1","K2",
             "e11","e12","e21","e22",
             "h11","h12","h21","h22",
             "a11","a12","a21","a22",
             "m1","m2","w11","w22")

# The probability quantile function used to draw values of each parameter. Since I did not have information on the distribution of parameters I used a uniform distribution to have equiprobable values in a range from min to max:
q <- replicate(24, "qunif")

# A list of the parameters of the above probability distributions. For the uniform distribution the parameters are minimum and maximum values. So in this case each element of the list is another list with the max and min of the corresponding model parameter, in the same order as above.
# for 4 species (R1,R2,C1,C2)
q.arg.4 <- list(list(min= 0.5, max = 2),list(min = 0.5, max = 2),list(min = 0.5, max = 2),list(min = 0.5, max = 2),
              list(min = 0.5, max = 2),list(min = 0.5, max = 2),list(min = 2, max = 5),list(min = 2, max = 5),
              list(min = 0.5, max = 0.9),list(min = 0.5, max = 0.9),list(min = 0.5, max = 0.9),list(min = 0.5, max = 0.9),
              list(min = 0.1, max = 0.6),list(min = 0.1, max = 0.6),list(min = 0.1, max = 0.6),list(min = 0.1, max = 0.6),
              list(min = 0.1, max = 4),list(min = 0.1, max = 4),list(min = 0.1, max = 4),list(min = 0.1, max = 4),
              list(min = 0.1, max = 1),list(min = 0.1, max = 1),list(min = 0.05, max = 0.95),list(min = 0.05, max = 0.95))

## Run simulation ----

## 4. Now create the hypercube and run the model for each combination of parameters with LHS, the working horse function of pse package. This is for the 4-species model
system.time(myLHS.4 <- LHS(model=modelRun, factors=factors, N=n.p, q=q, q.arg=q.arg.4, nboot=0)) # I used nboot=0 because I did not need significance of partial correlations

## Confidence bars crossing the zero horizontal line indicate no-significant partial correlation
#plotprcc(myLHS.4, ylab = c("R1","R2","C1","C2","steady")) 

## 5. Get the table of parameter combinations with get.data
hypercube.4 <- get.data(myLHS.4)

## 6. and the output for each combination with get.results
hypercube.4 <- cbind(hypercube.4, get.results(myLHS.4)) # appending columns of outputs to the hypercube matrix
names(hypercube.4)[25:29] <- c("R1", "R2", "C1", "C2", "steady") # cosmetic

## save data
write.csv(x = hypercube.4, file = "data/hypercube.4sp.csv")

