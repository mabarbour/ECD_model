
## Description: Ensure feasibility and stability of parameters for both 4 species and 3 species models

## source general parameters and functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

## load data ----
df.4sp.feas <- read.csv("data/hypercube.4sp.csv") %>% 
  tbl_df() %>% 
  select(-X) %>%
  filter(R1 > abund.thres, R2 > abund.thres, C1 > abund.thres, C2 > abund.thres, steady == 1)
df.4sp.feas # 525 feasible equilibriums where all 4 species coexist

## double-check stability and dynamics of 4 species equilibrium (negative real eigenvalue and presence of imaginary eigenvalue) ----
stable.4sp <- c()
spin.4sp <- c()
for(i in 1:dim(df.4sp.feas)[1]){
  set <- as.matrix(df.4sp.feas[i,5:28])
  eq.jac.set <- jacobian.full(y = set[1,c("R1","R2","C1","C2")], func = ECD_model, parms = set[1,1:20])
  stable.4sp[i] <- max(Re(eigen(eq.jac.set)$values))
  spin.4sp[i] <- max(Im(eigen(eq.jac.set)$values))
}

df.4sp.feas$max.Re.eigen <- stable.4sp
df.4sp.feas$max.Im.eigen <- spin.4sp

length(which(spin.4sp == 0)) # 8 of them appear to have monotonic dynamics.
length(which(stable.4sp < 0)) # 440 feasible and stable equilibriums with 4 sp. coexisting

df.4sp.feas.stab <- filter(df.4sp.feas, max.Re.eigen < 0)

## determine which parameter sets also have feasible equilibriums if one consumer species is removed ----
df <- df.4sp.feas.stab
feas.3sp.list <- list()
for(i in 1:dim(df)[1]){
  
  ## specify the simulation data to use
  set <- as.matrix(df[i, ]) 
  
  ## use initial state values
  C1.states <- set[1,c("R1_0","R2_0","C1_0")]
  names(C1.states) <- c("R1","R2","C1")
  C2.states <- set[1,c("R1_0","R2_0","C2_0")]
  names(C2.states) <- c("R1","R2","C2")
  
  ## grab parameter sets
  C1.pset <- set[1,c("r1","r2","K1","K2","e11","e12","h11","h12","a11","a12","m1","w11")]
  C2.pset <- set[1,c("r1","r2","K1","K2","e21","e22","h21","h22","a21","a22","m2","w22")]
  
  ## calculate new equilibrium values for 3 species
  C1.out <- safe.runsteady(y = C1.states, func = ECD_model.C1.3sp, parms = C1.pset)
  C2.out <- safe.runsteady(y = C2.states, func = ECD_model.C2.3sp, parms = C2.pset)
  
  ## calculate Jacobian matrix with new equilibrium values  
  
  # C1
  if(any(is.na(C1.out$y)) == FALSE){
    C1.eq.jac <- jacobian.full(y = C1.out$y, func = ECD_model.C1.3sp, parms = C1.pset) 
    C1.max.eigen <- max(Re(eigen(C1.eq.jac)$values))
  }
  if(any(is.na(C1.out$y)) == TRUE){
    C1.max.eigen <- NA
  }
  
  # C2
  if(any(is.na(C2.out$y)) == FALSE){
    C2.eq.jac <- jacobian.full(y = C2.out$y, func = ECD_model.C2.3sp, parms = C2.pset)
    C2.max.eigen <- max(Re(eigen(C2.eq.jac)$values))
  }
  if(any(is.na(C2.out$y)) == TRUE){
    C2.max.eigen <- NA
  }
  
  ## output all data, including stability (maximum real eigenvalue) 
  feas.3sp.list[[i]] <- c(sim.number = i,
                          C1.R1 = as.numeric(C1.out$y["R1"]),
                          C1.R2 = as.numeric(C1.out$y["R2"]),
                          C1.C1 = as.numeric(C1.out$y["C1"]),
                          C2.R1 = as.numeric(C2.out$y["R1"]),
                          C2.R2 = as.numeric(C2.out$y["R2"]),
                          C2.C2 = as.numeric(C2.out$y["C2"]),
                          steady.C1.3sp = attr(C1.out, "steady"),
                          steady.C2.3sp = attr(C2.out, "steady"),
                          C1.max.eigen = C1.max.eigen,
                          C2.max.eigen = C2.max.eigen)
}

all.feas.stable.df <- cbind.data.frame(df.4sp.feas.stab, ldply(feas.3sp.list)) %>%
  filter(C1.R1 > abund.thres, C1.R2 > abund.thres, C1.C1 > abund.thres, 
         C2.R1 > abund.thres, C2.R2 > abund.thres, C2.C2 > abund.thres,
         steady.C1.3sp == 1, steady.C2.3sp == 1,
         C1.max.eigen < 0, C2.max.eigen < 0)

## Write out the simulation data
write.csv(all.feas.stable.df, file = "data/all.feas.stable.df.csv")

## Check certain parameter sets ----
row <- 222
ck.set <- as.matrix(df.4sp.feas.stab[row, ])

## 4 sp. scenario
ck.4.states <- ck.set[1,c("R1_0","R2_0","C1_0","C2_0")]
names(ck.4.states) <- c("R1","R2","C1","C2")
ck.4.pset <- ck.set[1,5:24]

safe.runsteady(y = ck.4.states, func = ECD_model, parms = ck.4.pset, stol = 1e-5)
matplot.0D(ode(ck.4.states, 1:1000, ECD_model, ck.4.pset))

## 3 sp. scenarios
ck.C1.states <- ck.set[1,c("R1_0","R2_0","C1_0")]
names(ck.C1.states) <- c("R1","R2","C1")
ck.C2.states <- ck.set[1,c("R1_0","R2_0","C2_0")]
names(ck.C2.states) <- c("R1","R2","C2")
ck.C1.pset <- ck.set[1,c("r1","r2","K1","K2","e11","e12","h11","h12","a11","a12","m1","w11")]
ck.C2.pset <- ck.set[1,c("r1","r2","K1","K2","e21","e22","h21","h22","a21","a22","m2","w22")]

# C1
safe.runsteady(y = ck.C1.states, func = ECD_model.C1.3sp, parms = ck.C1.pset)
matplot.0D(ode(ck.C1.states, 1:1000, ECD_model.C1.3sp, ck.C1.pset))
jacobian.full(y = C1.states, func = ECD_model.C1.3sp, parms = C1.pset)

# C2
safe.runsteady(y = ck.set[1,c("R1","R2","C2")], func = ECD_model.C2.3sp, parms = ck.C2.pset)
matplot.0D(ode(ck.C2.states, 1:1000, ECD_model.C2.3sp, ck.C2.pset))


