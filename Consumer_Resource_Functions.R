

#### CONSUMER-RESOURCE MODELS ####

# note that all functions must be in the necessary format "function(Time,State,Pars)" for solving with the ode() function in the R package "deSolve"

#### MacArthur 1972, Geographical Ecology ####

## 4 species model
MacArthur_2C_2R <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)), {
    
    ## Consumer functional responses
    C1R1fxn <- a11 * R1 
    C1R2fxn <- a12 * R2 
    C2R1fxn <- a21 * R1 
    C2R2fxn <- a22 * R2 
    
    ## Continuous-time dynamics 
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - C2 * C2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - C2 * C2R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## 3 species model: C1, R1, R2
MacArthur_1C_2R <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)), {
    
    ## Consumer functional responses
    C1R1fxn <- a11 * R1 
    C1R2fxn <- a12 * R2 
    
    ## Continuous-time dynamics
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn 
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn 
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt))) 
  })
}


#### Lawlor and Smith 1976, Am Nat ####

# Extension of MacArthur consumer-resource model with habitat preference.
# Assumes that resources occur in distinct habitats and the consumers spend a proportion of time 'w' in each habitat foraging for the resource. 
# Note that wij = (1 - wii), but I've left it general in the equation

## 4 species model
LawlorSmith_2C_2R <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)), {
    
    ## Consumer functional responses
    C1R1fxn <- w11 * a11 * R1 
    C1R2fxn <- w12 * a12 * R2 
    C2R1fxn <- w21 * a21 * R1 
    C2R2fxn <- w22 * a22 * R2 
    
    ## Continuous-time dynamics 
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - C2 * C2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - C2 * C2R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## 3 species model: C1, R1, R2
LawlorSmith_1C_2R <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)), {
    
    ## Consumer functional responses
    C1R1fxn <- w11 * a11 * R1 
    C1R2fxn <- w12 * a12 * R2 
    
    ## Continuous-time dynamics
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn 
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn 
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt))) 
  })
}

#### McCann 2005 ####

## 4 species model: 2 consumers, 2 resources
McCann_2C_2R <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)), {
      
    ## Consumer preference functions
    W11 <- (w11 * R1)/(w11 * R1 + w12 * R2)
    W12 <- (w12 * R2)/(w11 * R1 + w12 * R2)
    W21 <- (w21 * R1)/(w21 * R1 + w22 * R2)
    W22 <- (w22 * R2)/(w21 * R1 + w22 * R2)

    ## Consumer functional responses
    C1R1fxn <- (W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- (W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C2R1fxn <- (W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    C2R2fxn <- (W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    ## Continuous-time dynamics
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - C2 * C2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - C2 * C2R2fxn
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## 3 species model: C1, R1, R2
McCann_1C_2R <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    ## Consumer preference functions
    W11 <- (w11 * R1)/(w11 * R1 + w12 * R2)
    W12 <- (w12 * R2)/(w11 * R1 + w12 * R2)
    
    ## Consumer functional responses
    C1R1fxn <- (W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- (W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    
    ## Continuous-time dynamics
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn 
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn 
    dC1.dt <- C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt)))
  })
}
