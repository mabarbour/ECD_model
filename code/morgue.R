# stode() identifies the steady-state of the model. It appears to be a lot of faster, and actually more accurate than runsteady...
safe.stode <- failwith(default = c(R1 = NA, R2 = NA, C1 = NA, C2 = NA, steady = 0), stode)

## using a threshold value to identify equilibrium values (roots) ----
rootfun <- function(Time, State, Pars){
  dstate <- unlist(ECD_model(Time, State, Pars)) # rate of change
  return(sum(abs(dstate)) - 1e-8) # threshold difference
}

### MAYBE USEFUL
set <- as.matrix(df[14, ])
C1.pset <- set[1,c("r1","r2","K1","K2","e11","e12","h11","h12","a11","a12","m1","w11")]
C2.pset <- set[1,c("r1","r2","K1","K2","e21","e22","h21","h22","a21","a22","m2","w22")]
dynamic_matplot(param.vector = C2.pset, 
                init.state = set[1,c("R1","R2","C2")], 
                model = ECD_model.C2.3sp, sim.length = 1000, ylim = c(0,5))

stode(y = set[1,c("R1","R2","C2")],func = ECD_model.C2.3sp, parms = C2.pset, positive = TRUE)

#### SHOULDN'T NEED THIS BELOW UNLESS I WANT TO SEE IF 4 SP AND 3 SP HAVE SIMILAR FEASIBILITY/STABILITY SPACE.
## Repeat the same thing but for 3 species ----


## Setup simulation ----

## Create a function that gets the model parameters and time sequence and outputs the final size of populations
oneRun.3sp <- function(R1_0, R2_0, C1_0, #C2_0,
                       r1, r2, K1, K2,
                       e11, e12, #e21, e22,
                       h11, h12, #h21, h22,
                       a11, a12, #a21, a22,
                       m1, w11){  # m2, w22
  # give parameters, initial state variables, and simulation length
  parameters <- c(r1 = r1, r2 = r2, K1 = K1, K2 = K2,
                  e11 = e11, e12 = e12, #e21 = e21, e22 = e22,
                  h11 = h11, h12 = h12, #h21 = h21, h22 = h22,
                  a11 = a11, a12 = a12, #a21 = a21, a22 = a22,
                  m1 = m1, w11 = w11) #, m2 = m2, w22 = w22
  state <- c(R1 = R1_0, R2 = R2_0, C1 = C1_0) #, C2 = C2_0
  times <- c(0, Inf) # length of simulation
  
  # using failwith() to make function give output even with a fatal error
  # runsteady() identifies the steady-state of the model.
  #safe.runsteady <- failwith(default = c(R1 = NA, R2 = NA, C1 = NA, C2 = NA, steady = 0), runsteady)
  
  out <- safe.stode(y = state, func = ECD_model.3sp, parms = parameters, positive = TRUE) #times = times, 
  
  # return steady states
  if(class(out) == "list"){
    return(c(out$y, steady = attr(out, "steady")))
  }
  if(class(out) == "numeric"){
    return(out)
  }
}

## Create a function with mapply that gets a matrix of parameter combinations (combinations are lines and parameters are columns) and returns the output of the function above for each parameter combination.
modelRun.3sp <- function (my.pars) {
  return(mapply(oneRun.3sp, 
                my.pars[,1], my.pars[,2], my.pars[,3], #my.pars[,4], 
                my.pars[,4], my.pars[,5], my.pars[,6], my.pars[,7],
                my.pars[,8], my.pars[,9], #my.pars[,11], my.pars[,12],
                my.pars[,10], my.pars[,11], #my.pars[,15], my.pars[,16],
                my.pars[,12], my.pars[,13], #my.pars[,19], my.pars[,20],
                my.pars[,14], my.pars[,15])) #my.pars[,22],my.pars[,24]
}

## Now prepare the hypercube.

# You have to specify a character vector to name the state variables and the model parameters
factors.3sp <- c("R1_0","R2_0","C1_0",  #"C2_0",
                 "r1","r2","K1","K2",
                 "e11","e12",   #"e21","e22",
                 "h11","h12",   #"h21","h22",
                 "a11","a12",   #"a21","a22",
                 "m1","w11")   #"m2","w22"

# The probability quantile function used to draw values of each parameter. Since I did not have information on the distribution of parameters I used a uniform distribution to have equiprobable values in a range from min to max:
q.3sp <- replicate(15, "qunif")
# Specificy the number of parameter sets to generate
n.p.3sp <- 500

# parameter list for 3 species (R1,R2,C1). Set C2 = 0 and relevant parameters for C2 = 0.
q.arg.3 <- list(list(min= 0.5, max = 2),list(min = 0.5, max = 2),list(min = 0.5, max = 2), #list(min = 0, max = 0),
                list(min = 0.5, max = 2),list(min = 0.5, max = 2),list(min = 2, max = 5),list(min = 2, max = 5),
                list(min = 0.5, max = 0.9),list(min = 0.5, max = 0.9),#list(min = 0, max = 0),list(min = 0, max = 0),
                list(min = 0.1, max = 0.6),list(min = 0.1, max = 0.6),#list(min = 0, max = 0),list(min = 0, max = 0),
                list(min = 0.1, max = 4),list(min = 0.1, max = 4),#list(min = 0, max = 0),list(min = 0, max = 0),
                list(min = 0.1, max = 1),list(min = 0.05, max = 0.95)) #list(min = 0.1, max = 0) ,list(min = 0, max = 0)


## Create hypercube and run model for 3-species (R1,R2,C1) 
myLHS.3 <- LHS(model=modelRun.3sp, factors=factors.3sp, N=n.p.3sp, q=q.3sp, q.arg=q.arg.3, nboot=0) 

## Get the table of parameter combinations with get.data
hypercube.3 <- get.data(myLHS.3)

## and the output for each combination with get.results
hypercube.3 <- cbind(hypercube.3, get.results(myLHS.3)) # appending columns of outputs to the hypercube matrix
names(hypercube.3)[16:19] <- c("R1", "R2", "C1", "steady") # cosmetic #"C2"

## save data
write.csv(x = hypercube.3, file = "Google Drive/ms_ecd/hypercube.3sp.stode.csv")


######## everything below may be useful but not for getting the simulation data
pos.stable.4 <- filter(hypercube.4, R1 > 1e-3, R2 > 1e-3, C1 > 1e-3, C2 > 1e-3, steady == 1)
row <- 6
dynamic_matplot(param.vector = pos.stable.4[row,5:24], init.state = c(R1 = pos.stable.4[row,"R1_0"],
                                                                      R2 = pos.stable.4[row,"R2_0"],
                                                                      C1 = pos.stable.4[row,"C1_0"],
                                                                      C2 = pos.stable.4[row,"C2_0"]), sim.length = 500, model = ECD_model, ylim = c(0,5))
## 5. And now use exploratory statistics to identify interesting regions of the parameter space
## pse package provides some options like
## 5a. Cumulative distribution of outputs
#plotecdf(myLHS, stack=TRUE) # about 40% of the population sizes are ~zero, so coexistence is not granted
## 5b. Scatterplots of response x parameter values
#plotscatter(myLHS, add.lm=FALSE)
## Again we see that many runs resulted in exclusion of one of the species
## We see also that initial conditions do not affect the final result, but
## competition coefficients do.
## Partial rank correlations (run LHS with nboot to get confidence intervals)
## Which are the non-parametric correlation of each parameter with each output,
## with the effects of the other variables partialed out
## Confidence bars crossing the zero horizontal line indicate no-significant partial correlation
#plotprcc(myLHS)
## Final population sizes have a strong correlation with competition coefficients




####### haven't explored below this
##7. And now you are free to use any exploratory tool to identify interesting regions of the parameter space
## There is not a single recipe, but a first guess is to try univariate and then bivariate and multivariate analyses.
## In some case you can also use linear models such as multiple regressions to identify patterns.

## We already know that there is coexistence in this parameter space.
## but in how many runs?
sum(hypercube$X>1e-6&hypercube$Y>1e-6)

## Which parameter combinations ensue coexistence?
## Univariate exploration with boxplot
## Create a logical variable to flag coexistence
hypercube$coexistence <- hypercube$X>1e-6&hypercube$Y>1e-6
## boxplots
par(mfrow=c(2,2))
boxplot(X0~coexistence, data=hypercube, xlab="Coexistence")
boxplot(Y0~coexistence, data=hypercube, xlab="Coexistence")
boxplot(a~coexistence, data=hypercube, xlab="Coexistence")
boxplot(b~coexistence, data=hypercube, xlab="Coexistence")
par(mfrow=c(1,1))
## There is an upper bound of coefficients to have coexistence

## Which combination of coefficients ensue coexistence?
## A bivariate scatterplot
plot(a~b, data=hypercube, type="n") ## to scale the plot for the whole parameter space
points(a~b, data=hypercube, subset=hypercube$X>1e-6&hypercube$Y>1e-6, col="blue")
points(a~b, data=hypercube, subset=!(hypercube$X>1e-6&hypercube$Y>1e-6), col="red")
## Coexistence possible only if both  coeficients < 1, roughly
## Further exploration can  reveal more!


### A more complicated example: Disease spread in a Suscetible-Infected-Recovered epidemic model ###

## 1. Create a function that receive model parameters and time sequence
## and outputs a logical value: disease has spread or not
oneRun.sir <- function(S0, I0, R0, r, a, time=seq(0, 50, by = 0.01)){
  ## parameters: transmission and recovering rates
  parameters <- c(r = r, a = a)
  ## initial population sizes: start with a single infected individual to test invasibility
  state <- c(S = S0, I = I0, R = R0)
  ## The function to be integrated
  sir <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      dS = -r*S*I
      dI = r*S*I - a*I
      dR = a*I
      list(c(dS, dI, dR))
    })
  }
  ## Integrating
  out <- ode(y = state, times = time, func = sir, parms = parameters)
  return(any(out[,3]>I0)) # returns a logical: any infected number larger than initial number?
}

## 2. Create a function that gets matrix of parameter combinations
## (combinations are lines and parameters are columns)
## and returns the output of the integration for each parameter combination.
## To do this, use mapply
## As we are interested in disease spread we restrict I0 to one and R0 to zero
modelRun.sir <- function (my.pars) {
  return(mapply(oneRun.sir, my.pars[,1], 1, 0, my.pars[,2], my.pars[,3]))
}

## 3. Now prepare the hypercube. You have to specify:
## A string vector to name the parameters
factors <- c("S0", "r", "a")
## The probability quantile function used to draw values of each parameter
## for equiprobable values use uniform distribution:
q <- c("qunif", "qunif", "qunif")
## A list of the parameters of the above probability distributions
## for the uniform the parameters are minimum and maximun values
## So in this case Each element of the list are is a list
## with the max and min of the corresponding model parameter
q.arg <- list(list(min=10, max=200), list(min=0.001, max=0.01), list(min=0, max=1))

## 4. Now create the hypercube and run the model for each combination of parameters
## 200 combinations of parameters
myLHS <- LHS(model=modelRun.sir, factors=factors, N=200, q=q, q.arg=q.arg)

## 5. you get the table of parameter combinations with get.data
hypercube <- get.data(myLHS)

## 6. and the output for each combination with get.results
hypercube$Spread <- get.results(myLHS)

## Now explore the relationships between parameters and outputs
## Start looking at the effect of each parameter
## In this case a boxplot of parameter values by output
par(mfrow=c(1,3))
boxplot(S0~Spread, data=hypercube, xlab="Disease spread", ylab="Initial number of susceptibles")
boxplot(r~Spread, data=hypercube, xlab="Disease spread", ylab="Transmission rate")
boxplot(a~Spread, data=hypercube, xlab="Disease spread", ylab="Recovering rate")
par(mfrow=c(1,1))

## You can also explora relatioships between parameters conditioned
## conditioned to the output
## Relationship a X r when disease spread or not
## only populations with initial size < median of all initial size
plot(a~r, data=hypercube, type="n") # to show all ranges of the parameters
points(a~r, data=hypercube, subset=S0<median(S0)&Spread==1, col="red")
points(a~r, data=hypercube, subset=S0<median(S0)&Spread==0, col="blue")

## And you can also explore relationships between ratios of two parameters
## in function of additional  parameters
plot(I(a/r)~S0, data=hypercube, type="n") # to show all ranges of the parameters
points(I(a/r)~S0, subset=Spread==1, data=hypercube, col="red")
points(I(a/r)~S0, subset=Spread==0, data=hypercube, col="blue")
abline(0,1) ## equivalence line
## BINGO! Disease spreads if S0 > a/r

#######  NEED TO ALTER 3 SP SIMS SO I CAN EVENTUALLY CALCULATE DIVERGENCE.

## C1 3 species: Setup For-Loop for simulation ----
mut.reps <- 500 # number of mutation iterations
sim.list.C1 <- list()
for(i in 1:dim(df.4sp.feas)[1]){
  #i <- 1
  # choose parameter set to examine and calculate the maximum real eigenvalue
  set <- as.matrix(df.4sp.feas[i, 5:29])#[ ,-c(7:8,11:12,15:16,18,20,24)]
  #set <- set[1, -c(7:8,11:12,15:16,18,20,24)] #<- 0
  set[1, c("e21","e22","h21","h22","a21","a22","m2","w22","C2")] <- 0
  
  ## calculate new equilibrium values for 3 species
  out <- safe.stode(y = set[1,c("R1","R2","C1")], func = ECD_model.C1.3sp, parms = set[1,c(1:6,9:10,13:14,17,19)], positive = TRUE)
  
  sub.sim.list.C1 <- list()
  if(any(out$y < 1e-3) == FALSE){
    eq.jac.set <- jac.norm.calc(data = c(set[1,1:20],out$y,C2=0), jac.mat = jac)
    
    ## setup and run sub-simulations
    
    sub.sim.list.C1[[1]] <- c(set[1,1:20],out$y,C2=0, steady = attr(out, "steady"), max.eigen = max(Re(eigen(eq.jac.set[1:3,1:3])$values)), sp = NA, mut.suc = NA, p.try = NA) # input first data vector
  }
  if(any(out$y < 1e-3) == TRUE){
    sub.sim.list.C1[[1]] <- c(set[1,1:20],out$y,C2=0, steady = 0.5, max.eigen = NA, species = NA, mut.suc = NA, p.try = NA)
  }
  
  for(j in 2:mut.reps){
    
    #browser()
    #j <- 2
    # if the prior simulation is no longer at a steady state, just return NAs for the rest of the mutations reps
    if(sub.sim.list.C1[[j-1]][25] < 1){
      
      sub.sim.list.C1[[j]] <- c(sub.sim.list.C1[[j-1]][1:25], max.eigen = NA, species = NA, mut.suc = NA, p.try = NA)
      
    }
    
    # if the prior simulation is at a steady state
    # if all positive, proceed with simulation
    if(sub.sim.list.C1[[j-1]][25] == 1 & any(sub.sim.list.C1[[j-1]][21:23] < 1e-3) == FALSE){
      
      # extract parameters and state variables
      parameters <- sub.sim.list.C1[[j-1]][1:20]
      state <- sub.sim.list.C1[[j-1]][21:24]
      
      # assign mutation in attack rate to one of the species
      sp <- 1 # assign to C1
      mut <- ifelse(rbinom(1,1,0.5) == 0, -0.1, 0.1)
      
      # determine whether mutant can invade or not by calculating whether it has positive growth rate (i.e. positive eigen value)
      new.ps <- parameters
      #browser()
      if(sp == 1){
        new.ps["a11m"] <- new.ps["a11"] + mut
        new.ps["a12m"] <- new.ps["a12"] - mut
        new.ps["a21m"] <- 0
        new.ps["a22m"] <- 0
        jac.mut <- jac.mut.calc(data = c(new.ps, state), jac.mat = jac)
        eigen.mut <- jac.mut[5,5] # max(Re(eigen(jac.mut)$values))
      }
      if(sp == 2){
        new.ps["a11m"] <- 0
        new.ps["a12m"] <- 0
        new.ps["a22m"] <- new.ps["a22"] + mut
        new.ps["a21m"] <- new.ps["a21"] - mut
        jac.mut <- jac.mut.calc(data = c(new.ps, state), jac.mat = jac)
        eigen.mut <- jac.mut[6,6] # max(Re(eigen(jac.mut)$values))
      }
      
      #browser()
      # if the mutant can invade and attack rates are positive, update attack rates and print new equilibrium values
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
        #out <- safe.runsteady(y = state, times = c(0,Inf), func = ECD_model, parms = parameters)
        out <- safe.stode(y = state[1:3], func = ECD_model.C1.3sp, parms = parameters[-c(7:8,11:12,15:16,18,20)], positive = TRUE)
        
        # return steady states and maximum real eigenvalues
        if(class(out) == "list" & any(out$y < 1e-3) == FALSE){
          #if(class(out) == "list" & all(out$y > 0) == TRUE){
          
          eq.jac <- jac.norm.calc(data = c(parameters, out$y, C2 = 0), jac.mat = jac)
          
          sub.sim.list.C1[[j]] <- c(parameters, out$y, C2 = 0, steady = attr(out, "steady"), max.eigen = max(Re(eigen(eq.jac[1:3,1:3])$values)), species = sp, mut.suc = 1, p.try = 1)
          
        }
        if(class(out) == "list" & any(out$y < 1e-3) == TRUE){
          #if(class(out) == "list" & all(out$y > 0) == FALSE){
          
          eq.jac <- jac.norm.calc(data = c(parameters, out$y, C2 = 0), jac.mat = jac)
          
          sub.sim.list.C1[[j]] <- c(parameters, out$y, C2 = 0, steady = 0.5, max.eigen = NA, species = sp, mut.suc = 1, p.try = 1)
          
        }
        
        #if(class(out) == "numeric"){
        # sub.sim.list.4[[j]] <- c(parameters, out, max.eigen = NA, species = sp)
        #}
      }
      
      # if the mutant can't invade or the evolved attack rates are negative, keep attack rates the same and print the same equilibrium values with maximum real eigenvalue
      if(eigen.mut < 0 & all(new.ps >= 0) == TRUE){
        
        eq.jac <- jac.norm.calc(data = c(parameters, state), jac.mat = jac)
        
        sub.sim.list.C1[[j]] <- c(parameters, state, steady = sub.sim.list.C1[[j-1]][25], max.eigen = max(Re(eigen(eq.jac[1:3,1:3])$values)), species = sp, mut.suc = 0, p.try = 1)
        
      }
      if(all(new.ps >= 0) == FALSE){
        
        eq.jac <- jac.norm.calc(data = c(parameters, state), jac.mat = jac)
        
        sub.sim.list.C1[[j]] <- c(parameters, state, steady = sub.sim.list.C1[[j-1]][25], max.eigen = max(Re(eigen(eq.jac[1:3,1:3])$values)), species = sp, mut.suc = NA, p.try = -1)
        
      }
    }
  }
  # turn sub-simulation into data frame and output to simulation list.
  sim.list.C1[[i]] <- sub.sim.list.C1 %>% ldply()
}


## convert list into data frame, add sequence and simulation number, then print out data
evol.sim.df.C1 <- sim.list.C1 %>% ldply() %>% 
  mutate(sequence = rep(1:mut.reps, times = dim(df.4sp.feas)[1]), 
         sim.number = rep(1:dim(df.4sp.feas)[1], each = mut.reps)) 

write.csv(x = evol.sim.df.C1, file = "Google Drive/ms_ecd/evol.sim.df.C1.3sp.csv")

## C2 3 species: Setup For-Loop for simulation ----
# took a couple of hours to complete with 50 mutation reps, probably want to boost it to at least 1000.

mut.reps <- 500 # number of mutation iterations
sim.list.C2 <- list()
for(i in 1:dim(df.4sp.feas)[1]){
  #i <- 1
  # choose parameter set to examine and calculate the maximum real eigenvalue
  set <- as.matrix(df.4sp.feas[i, 5:29])#[ ,-c(7:8,11:12,15:16,18,20,24)]
  #set <- set[1, -c(7:8,11:12,15:16,18,20,24)] #<- 0
  set[1, c("e11","e12","h11","h12","a11","a12","m1","w11","C1")] <- 0
  
  ## calculate new equilibrium values for 3 species
  out <- safe.stode(y = set[1,c("R1","R2","C2")], func = ECD_model.C2.3sp, parms = set[1,c(1:4,7:8,11:12,15:16,18,20)], positive = TRUE)
  
  sub.sim.list.C2 <- list()
  if(any(out$y < 1e-3) == FALSE){
    eq.jac.set <- jac.norm.calc(data = c(set[1,1:20],out$y[1:2],C1=0,out$y[3]), jac.mat = jac)
    
    ## setup and run sub-simulations
    sub.sim.list.C2[[1]] <- c(set[1,1:20],out$y[1:2],C1=0,out$y[3], steady = attr(out, "steady"), max.eigen = max(Re(eigen(eq.jac.set[c(1,2,4),c(1,2,4)])$values)), sp = NA, mut.suc = NA, p.try = NA) # input first data vector
  }
  if(any(out$y < 1e-3) == TRUE){
    sub.sim.list.C2[[1]] <- c(set[1,1:20],out$y[1:2],C1=0,out$y[3], steady = 0.5, max.eigen = NA, species = NA, mut.suc = NA, p.try = NA)
  }
  #browser()
  for(j in 2:mut.reps){
    
    #browser()
    #j <- 2
    # if the prior simulation is no longer at a steady state, just return NAs for the rest of the mutations reps
    if(sub.sim.list.C2[[j-1]][25] < 1){
      
      sub.sim.list.C2[[j]] <- c(sub.sim.list.C2[[j-1]][1:25], max.eigen = NA, species = NA, mut.suc = NA, p.try = NA)
      
    }
    
    # if the prior simulation is at a steady state
    # if all positive, proceed with simulation
    if(sub.sim.list.C2[[j-1]][25] == 1 & any(sub.sim.list.C2[[j-1]][c(21,22,24)] < 1e-3) == FALSE){
      
      # extract parameters and state variables
      parameters <- sub.sim.list.C2[[j-1]][1:20]
      state <- sub.sim.list.C2[[j-1]][21:24]
      
      # assign mutation in attack rate to one of the species
      sp <- 2 # assign to C2
      mut <- ifelse(rbinom(1,1,0.5) == 0, -0.1, 0.1)
      
      # determine whether mutant can invade or not by calculating whether it has positive growth rate (i.e. positive eigen value)
      new.ps <- parameters
      #browser()
      if(sp == 1){
        new.ps["a11m"] <- new.ps["a11"] + mut
        new.ps["a12m"] <- new.ps["a12"] - mut
        new.ps["a21m"] <- 0
        new.ps["a22m"] <- 0
        jac.mut <- jac.mut.calc(data = c(new.ps, state), jac.mat = jac)
        eigen.mut <- jac.mut[5,5] # max(Re(eigen(jac.mut)$values))
      }
      if(sp == 2){
        new.ps["a11m"] <- 0
        new.ps["a12m"] <- 0
        new.ps["a22m"] <- new.ps["a22"] + mut
        new.ps["a21m"] <- new.ps["a21"] - mut
        jac.mut <- jac.mut.calc(data = c(new.ps, state), jac.mat = jac)
        eigen.mut <- jac.mut[6,6] # max(Re(eigen(jac.mut)$values))
      }
      
      #browser()
      # if the mutant can invade and attack rates are positive, update attack rates and print new equilibrium values
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
        #out <- safe.runsteady(y = state, times = c(0,Inf), func = ECD_model, parms = parameters)
        out <- safe.stode(y = state[c(1,2,4)], func = ECD_model.C2.3sp, parms = parameters[-c(5:6,9:10,13:14,17,19)], positive = TRUE)
        
        # return steady states and maximum real eigenvalues
        if(class(out) == "list" & any(out$y < 1e-3) == FALSE){
          #if(class(out) == "list" & all(out$y > 0) == TRUE){
          
          eq.jac <- jac.norm.calc(data = c(parameters, out$y[1:2], C1 = 0, out$y[3]), jac.mat = jac)
          
          sub.sim.list.C2[[j]] <- c(parameters, out$y[1:2], C1 = 0, out$y[3], steady = attr(out, "steady"), max.eigen = max(Re(eigen(eq.jac[c(1,2,4),c(1,2,4)])$values)), species = sp, mut.suc = 1, p.try = 1)
          
        }
        if(class(out) == "list" & any(out$y < 1e-3) == TRUE){
          #if(class(out) == "list" & all(out$y > 0) == FALSE){
          
          eq.jac <- jac.norm.calc(data = c(parameters, out$y[1:2], C1 = 0, out$y[3]), jac.mat = jac)
          
          sub.sim.list.C2[[j]] <- c(parameters, out$y[1:2], C1 = 0, out$y[3], steady = 0.5, max.eigen = NA, species = sp, mut.suc = 1, p.try = 1)
          
        }
        
        #if(class(out) == "numeric"){
        # sub.sim.list.4[[j]] <- c(parameters, out, max.eigen = NA, species = sp)
        #}
      }
      
      # if the mutant can't invade or the evolved attack rates are negative, keep attack rates the same and print the same equilibrium values with maximum real eigenvalue
      if(eigen.mut < 0 & all(new.ps >= 0) == TRUE){
        
        eq.jac <- jac.norm.calc(data = c(parameters, state), jac.mat = jac)
        
        sub.sim.list.C2[[j]] <- c(parameters, state, steady = sub.sim.list.C2[[j-1]][25], max.eigen = max(Re(eigen(eq.jac[c(1,2,4),c(1,2,4)])$values)), species = sp, mut.suc = 0, p.try = 1)
        
      }
      if(all(new.ps >= 0) == FALSE){
        
        eq.jac <- jac.norm.calc(data = c(parameters, state), jac.mat = jac)
        
        sub.sim.list.C2[[j]] <- c(parameters, state, steady = sub.sim.list.C2[[j-1]][25], max.eigen = max(Re(eigen(eq.jac[c(1,2,4),c(1,2,4)])$values)), species = sp, mut.suc = NA, p.try = -1)
        
      }
    }
  }
  # turn sub-simulation into data frame and output to simulation list.
  sim.list.C2[[i]] <- sub.sim.list.C2 %>% ldply()
}


## convert list into data frame, add sequence and simulation number, then print out data
evol.sim.df.C2 <- sim.list.C2 %>% ldply() %>% 
  mutate(sequence = rep(1:mut.reps, times = dim(df.4sp.feas)[1]), 
         sim.number = rep(1:dim(df.4sp.feas)[1], each = mut.reps)) 

write.csv(x = evol.sim.df.C2, file = "Google Drive/ms_ecd/evol.sim.df.C2.3sp.csv")

## maybe useful
## pick a parameter set and explore the dynamics
set <- as.matrix(df.4sp.unstable[50, ])[1, ]

dynamic_matplot(param.vector = set[5:24],
                init.state = c(R1 = as.numeric(set["R1_0"]), R2 = as.numeric(set[2]), C1 = as.numeric(set[3]), C2 = as.numeric(set[4])),
                sim.length = 1000,
                model = ECD_model, 
                ylim = c(0,max(set[c("K1","K2")])))
out.stode <- stode(y = c(R1 = as.numeric(set["R1_0"]), R2 = as.numeric(set[2]), C1 = as.numeric(set[3]), C2 = as.numeric(set[4])), func = ECD_model, parms = set[5:24], positive = TRUE)
out.run <- runsteady(y = c(R1 = as.numeric(set["R1_0"]), R2 = as.numeric(set[2]), C1 = as.numeric(set[3]), C2 = as.numeric(set[4])), func = ECD_model, parms = set[5:24])
