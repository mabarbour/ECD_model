######################################################
## Description: Dynamical model and other useful functions for checking dynamics
## Code author(s): Matt Barbour
## Email: barbour@zoology.ubc.ca
######################################################

## 4 species model: 2 consumers, 2 resources
Pachepsky_model <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    # constraints on consumer preference
    w12 <- (1 - w11)
    w21 <- (1 - w22)
      
    # consumer preference functions
    W11 <- (w11 * R1)/(w11 * R1 + w12 * R2)
    W12 <- (w12 * R2)/(w11 * R1 + w12 * R2)
    W21 <- (w21 * R1)/(w21 * R1 + w22 * R2)
    W22 <- (w22 * R2)/(w21 * R1 + w22 * R2)

    
    # consumer functional responses
    C1R1fxn <- (W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- (W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C2R1fxn <- (W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    C2R2fxn <- (W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn - C2 * C2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn - C2 * C2R2fxn
    dC1.dt <- -m1 * C1 # C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    dC2.dt <- -m2 * C2 # C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## 3 species model: C1, R1, R2
Pachepsky_model.C1.3sp <- function(Time,State,Pars) {
  # t,y,p is the necessary format for solving using the ode() function in R
  with(as.list(c(State,Pars)), {
    
    # constraints on consumer preference
    w12 <- (1 - w11)
    #w21 <- (1 - w22)
      
    # consumer preference functions
    W11 <- (w11 * R1)/(w11 * R1 + w12 * R2)
    W12 <- (w12 * R2)/(w11 * R1 + w12 * R2)
    #W21 <- (w21 * R1)/(w21 * R1 + w22 * R2)
    #W22 <- (w22 * R2)/(w21 * R1 + w22 * R2)
    
    
    # consumer functional responses
    C1R1fxn <- (W11 * a11 * R1)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    C1R2fxn <- (W12 * a12 * R2)/(1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
    #C2R1fxn <- (W21 * a21 * R1)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    #C2R2fxn <- (W22 * a22 * R2)/(1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
    
    # dynamical equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - C1 * C1R1fxn #- C2 * C2R1fxn
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - C1 * C1R2fxn #- C2 * C2R2fxn
    dC1.dt <- -m1 * C1 #C1 * (e11 * C1R1fxn + e12 * C1R2fxn - m1)
    #dC2.dt <- C2 * (e21 * C2R1fxn + e22 * C2R2fxn - m2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt))) #, dC2.dt
  })
}
