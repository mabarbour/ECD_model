
## Calculate Evolutionary Stable Strategies (ESS)
# Evolutionary optimum, where a mutant is no longer able to invade (first derivative = 0 and second derivative is negative)

## Load required libraries
library(rootSolve)
library(tidyverse)

## Source in required functions
source('identify_steady_state.R')
source('eco_evo_sim_function.R')
source('Consumer_Resource_Functions.R')
source('CR_dynamics.R')

## Custom functions

# calculate attack rate on alternative resource
get_a12 <- function(x, A, n) A*(1-(x/A)^n)^(1/n)

## MacArthur model
C1_ESS_2C_2R_MacArthur <- function(x, A, n, m, e, dx){
  
  # calculate a12 in terms of a11 (x; tradeoff)
  a11 = x
  a12 = A*(1-(a11/A)^n)^(1/n)
  
  # calculate Resource abundances at equilibrium
  eq.R1 = 1/(a11 + a12) * m/e
  eq.R2 = 1/(a11 + a12) * m/e
  
  # calculate mutant attack rates
  a11m = x + dx
  a12m = A*(1-(a11m/A)^n)^(1/n)
  
  # calculate mutant growth rate of C1 when rare
  e*a11m*eq.R1 + e*a12m*eq.R2 - m
}

C1_ESS_1C_2R_MacArthur <- function(x, A, n, m, e, dx, r, K, R1, R2, C1){
  
  # calculate a12 in terms of a11 (x; tradeoff)
  a11 <- x
  a12 <- A*(1-(a11/A)^n)^(1/n)
  
  # calculate Resource abundances at equilibrium
  get_steady <- safe.runsteady(y = c(R1 = R1, R2 = R2, C1 = C1), 
                               func = MacArthur_1C_2R, 
                               parms = c(a11 = a11, a12 = a12, 
                                         r1 = r, r2 = r,
                                         K1 = K, K2 = K,
                                         e11 = e, e12 = e,
                                         m1 = m)) 
  
  eq.R1 = get_steady$y["R1"][[1]]
  eq.R2 = get_steady$y["R2"][[1]]
  
  # calculate mutant attack rates
  a11m = x + dx
  a12m = A*(1-(a11m/A)^n)^(1/n)
  
  # calculate mutant growth rate of C1 when rare
  e*a11m*eq.R1 + e*a12m*eq.R2 - m
}

## General parameters
A <- 2
m <- 1
e <- 0.8

#### MacArthur, two consumer food webs

## Concave up tradeoff
n <- 0.85

# Plot
Mac_CU_neg <- curve(C1_ESS_2C_2R_MacArthur(x, A = A, n = n, m = m, e = e, dx = -1e-10), 0 + 0.01, A - 0.01)
Mac_CU_pos <- curve(C1_ESS_2C_2R_MacArthur(x, A = A, n = n, m = m, e = e, dx = 1e-10), 0 + 0.01, A - 0.01)
Mac_CU_df <- data.frame(x = c(Mac_CU_neg$x, Mac_CU_pos$x), 
                      y = c(Mac_CU_neg$y, Mac_CU_pos$y), 
                      dx = c(rep(-1e-10, length(Mac_CU_neg$x)), rep(1e-10, length(Mac_CU_pos$x))))

Mac_CU_root <- uniroot(C1_ESS_2C_2R_MacArthur, A = A, n = n, m = m, e = e, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
Mac_CU_root / (Mac_CU_root + get_a12(x = Mac_CU_root, A = A, n = n)) # mutant growth rate = 0 as complete generalist

# Indicative of complete specialization, since deviations from a complete generalist result in positive mutant growth rates
Mac_CU_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Mac_CU_root, linetype = "dotted") +
  ylab("Mutant growth rate")

# We can also see this by looking at how small deviations from complete generalist phenotype increase mutant growth rate
C1_ESS_2C_2R_MacArthur(x=Mac_CU_root, A = A, n = n, m = m, e = e, dx = 0.01)
C1_ESS_2C_2R_MacArthur(x=Mac_CU_root, A = A, n = n, m = m, e = e, dx = -0.01)


## Concave down tradeoff

n <- 1.15

# Plot
Mac_CD_neg <- curve(C1_ESS_2C_2R_MacArthur(x, A = A, n = n, m = m, e = e, dx = -1e-10), 0 + 0.01, A - 0.01)
Mac_CD_pos <- curve(C1_ESS_2C_2R_MacArthur(x, A = A, n = n, m = m, e = e, dx = 1e-10), 0 + 0.01, A - 0.01)
Mac_CD_df <- data.frame(x = c(Mac_CD_neg$x, Mac_CD_pos$x), 
                        y = c(Mac_CD_neg$y, Mac_CD_pos$y), 
                        dx = c(rep(-1e-10, length(Mac_CD_neg$x)), rep(1e-10, length(Mac_CD_pos$x))))

Mac_CD_root <- uniroot(C1_ESS_2C_2R_MacArthur, A = A, n = n, m = m, e = e, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
Mac_CD_root / (Mac_CD_root + get_a12(x = Mac_CD_root, A = A, n = n)) # mutant growth rate = 0 as complete generalist

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
Mac_CD_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Mac_CD_root, linetype = "dotted") +
  ylab("Mutant growth rate")

# We can also see this by looking at how small deviations from complete generalist phenotype decrease mutant growth rate
C1_ESS_2C_2R_MacArthur(x=Mac_CD_root, A = A, n = n, m = m, e = e, dx = 0.01)
C1_ESS_2C_2R_MacArthur(x=Mac_CD_root, A = A, n = n, m = m, e = e, dx = -0.01)


## Linear tradeoff

n <- 1

# Plot
Mac_LN_neg <- curve(C1_ESS_2C_2R_MacArthur(x, A = A, n = n, m = m, e = e, dx = -1e-10), 0 + 0.01, A - 0.01)
Mac_LN_pos <- curve(C1_ESS_2C_2R_MacArthur(x, A = A, n = n, m = m, e = e, dx = 1e-10), 0 + 0.01, A - 0.01)
Mac_LN_df <- data.frame(x = c(Mac_LN_neg$x, Mac_LN_pos$x), 
                        y = c(Mac_LN_neg$y, Mac_LN_pos$y), 
                        dx = c(rep(-1e-10, length(Mac_LN_neg$x)), rep(1e-10, length(Mac_LN_pos$x))))

Mac_LN_root <- uniroot(C1_ESS_2C_2R_MacArthur, A = A, n = n, m = m, e = e, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
Mac_LN_root / (Mac_LN_root + get_a12(x = Mac_LN_root, A = A, n = n)) # mutant growth rate = 0 as complete generalist

uniroot(C1_ESS_2C_2R_MacArthur, A = A, n = n, m = m, e = e, dx = 1e-10, interval = c(0 + 0.01, A - 0.01)) # no precision in estimate...

# Running into issues of machine precision...
Mac_LN_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Mac_LN_root, linetype = "dotted") +
  ylab("Mutant growth rate")

# We can also see this by looking at how small deviations from root do not affect mutant growth rate
C1_ESS_2C_2R_MacArthur(x=Mac_LN_root, A = A, n = n, m = m, e = e, dx = 0.01)
C1_ESS_2C_2R_MacArthur(x=Mac_LN_root, A = A, n = n, m = m, e = e, dx = -0.01)

# This indicates that for a linear tradeoff, all phenotypes are equally viable, and the amount of divergence depends on initial conditions.

#### MacArthur, one consumer food webs

# note that I use larger values of dx and sequences in between, because things get jumpy if 
# since I have to dynamically get the equilbrium resource values.

# need additional parameters because it is difficult to get an analytical solution with only one consumer
r <- 1
K <- 4
R1 <- 2
R2 <- 2
C1 <- 1

## Concave up tradeoff
n <- 0.85

# Plot
Mac_CU_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_MacArthur(x, A = A, n = n, m = m, e = e, r, K, R1, R2, C1, dx = -0.01))
Mac_CU_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_MacArthur(x, A = A, n = n, m = m, e = e, r, K, R1, R2, C1, dx = 0.01))
Mac_CU_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                             y = c(Mac_CU_oneC_neg, Mac_CU_oneC_pos), 
                             dx = c(rep(-0.01, length(Mac_CU_oneC_neg)), rep(0.01, length(Mac_CU_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
Mac_CU_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Mac_CU_root, linetype = "dotted") + # root for two consumer predicts root for one consumer
  ylab("Mutant growth rate")

## Concave down tradeoff
n <- 1.15

# Plot
Mac_CD_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_MacArthur(x, A = A, n = n, m = m, e = e, r, K, R1, R2, C1, dx = -0.01))
Mac_CD_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_MacArthur(x, A = A, n = n, m = m, e = e, r, K, R1, R2, C1, dx = 0.01))
Mac_CD_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                             y = c(Mac_CD_oneC_neg, Mac_CD_oneC_pos), 
                             dx = c(rep(-0.01, length(Mac_CD_oneC_neg)), rep(0.01, length(Mac_CD_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
Mac_CD_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Mac_CD_root, linetype = "dotted") + # root for two consumer predicts root for one consumer
  ylab("Mutant growth rate")


## Linear tradeoff
n <- 1

# Plot
Mac_LN_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_MacArthur(x, A = A, n = n, m = m, e = e, r, K, R1, R2, C1, dx = -0.01))
Mac_LN_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_MacArthur(x, A = A, n = n, m = m, e = e, r, K, R1, R2, C1, dx = 0.01))
Mac_LN_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                             y = c(Mac_LN_oneC_neg, Mac_LN_oneC_pos), 
                             dx = c(rep(-0.01, length(Mac_LN_oneC_neg)), rep(0.01, length(Mac_LN_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
Mac_LN_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = A/2, linetype = "dotted") + # root for two consumer predicts root for one consumer
  ylab("Mutant growth rate")


###################################################################################################################
#####  Lawlor and Smith 
###################################################################################################################

## Adding habitat preference
w <- 0.6

## Functions

C1_ESS_2C_2R_LawlorSmith <- function(x, A, n, m, e, w, dx){
  
  # calculate a12 in terms of a11 (x; tradeoff)
  a11 = x
  a12 = A*(1-(a11/A)^n)^(1/n)
  
  # calculate Resource abundances at equilibrium
  eq.R1 = 1/(w*a11 + (1-w)*a12) * m/e
  eq.R2 = 1/(w*a11 + (1-w)*a12) * m/e
  
  # calculate mutant attack rates
  a11m = x + dx
  a12m = A*(1-(a11m/A)^n)^(1/n)
  
  # calculate mutant growth rate of C1 when rare
  e*w*a11m*eq.R1 + e*(1-w)*a12m*eq.R2 - m
}

C1_ESS_1C_2R_LawlorSmith <- function(x, A, n, m, e, w, dx, r, K, R1, R2, C1){
  
  # calculate a12 in terms of a11 (x; tradeoff)
  a11 <- x
  a12 <- A*(1-(a11/A)^n)^(1/n)
  
  # calculate Resource abundances at equilibrium
  get_steady <- safe.runsteady(y = c(R1 = R1, R2 = R2, C1 = C1), 
                               func = LawlorSmith_1C_2R, 
                               parms = c(a11 = a11, a12 = a12, 
                                         r1 = r, r2 = r,
                                         K1 = K, K2 = K,
                                         e11 = e, e12 = e,
                                         w11 = w, w12 = (1-w),
                                         m1 = m)) 
  
  eq.R1 = get_steady$y["R1"][[1]]
  eq.R2 = get_steady$y["R2"][[1]]
  
  # calculate mutant attack rates
  a11m = x + dx
  a12m = A*(1-(a11m/A)^n)^(1/n)
  
  # calculate mutant growth rate of C1 when rare
  e*w*a11m*eq.R1 + e*(1-w)*a12m*eq.R2 - m
}

eq.R1_1C_2R_LawlorSmith <- function(x, A, n, m, e, w, r, K, R1, R2, C1){
  
  # calculate a12 in terms of a11 (x; tradeoff)
  a11 <- x
  a12 <- A*(1-(a11/A)^n)^(1/n)
  
  # calculate Resource abundances at equilibrium
  get_steady <- safe.runsteady(y = c(R1 = R1, R2 = R2, C1 = C1), 
                               func = LawlorSmith_1C_2R, 
                               parms = c(a11 = a11, a12 = a12, 
                                         r1 = r, r2 = r,
                                         K1 = K, K2 = K,
                                         e11 = e, e12 = e,
                                         w11 = w, w12 = (1-w),
                                         m1 = m)) 
  
  get_steady$y["R1"][[1]]
}

eq.R2_1C_2R_LawlorSmith <- function(x, A, n, m, e, w, r, K, R1, R2, C1){
  
  # calculate a12 in terms of a11 (x; tradeoff)
  a11 <- x
  a12 <- A*(1-(a11/A)^n)^(1/n)
  
  # calculate Resource abundances at equilibrium
  get_steady <- safe.runsteady(y = c(R1 = R1, R2 = R2, C1 = C1), 
                               func = LawlorSmith_1C_2R, 
                               parms = c(a11 = a11, a12 = a12, 
                                         r1 = r, r2 = r,
                                         K1 = K, K2 = K,
                                         e11 = e, e12 = e,
                                         w11 = w, w12 = (1-w),
                                         m1 = m)) 
  
  get_steady$y["R2"][[1]]
}

## Two consumer food webs

## Concave up tradeoff
n <- 0.85

# Plot
LS_CU_neg <- curve(C1_ESS_2C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, dx = -1e-10), 0 + 0.01, A - 0.01)
LS_CU_pos <- curve(C1_ESS_2C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, dx = 1e-10), 0 + 0.01, A - 0.01)
LS_CU_df <- data.frame(x = c(LS_CU_neg$x, LS_CU_pos$x), 
                        y = c(LS_CU_neg$y, LS_CU_pos$y), 
                        dx = c(rep(-1e-10, length(LS_CU_neg$x)), rep(1e-10, length(LS_CU_pos$x))))

LS_CU_root <- uniroot(C1_ESS_2C_2R_LawlorSmith, A = A, n = n, m = m, e = e, w = w, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
LS_CU_root / (LS_CU_root + get_a12(x = LS_CU_root, A = A, n = n)) # mutant growth rate = 0 as complete generalist

## REWRITE ## # Indicative of complete specialization, since deviations from a complete generalist result in positive mutant growth rates
LS_CU_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = LS_CU_root, linetype = "dotted") +
  ylab("Mutant growth rate")

## REWRITE ## # We can also see this by looking at how small deviations from complete generalist phenotype increase mutant growth rate
C1_ESS_2C_2R_LawlorSmith(x=LS_CU_root, A = A, n = n, m = m, e = e, w = w, dx = 0.01)
C1_ESS_2C_2R_LawlorSmith(x=LS_CU_root, A = A, n = n, m = m, e = e, w = w, dx = -0.01)


## Concave down tradeoff

n <- 1.15

# Plot
LS_CD_neg <- curve(C1_ESS_2C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, dx = -1e-10), 0 + 0.01, A - 0.01)
LS_CD_pos <- curve(C1_ESS_2C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, dx = 1e-10), 0 + 0.01, A - 0.01)
LS_CD_df <- data.frame(x = c(LS_CD_neg$x, LS_CD_pos$x), 
                        y = c(LS_CD_neg$y, LS_CD_pos$y), 
                        dx = c(rep(-1e-10, length(LS_CD_neg$x)), rep(1e-10, length(LS_CD_pos$x))))

LS_CD_root <- uniroot(C1_ESS_2C_2R_LawlorSmith, A = A, n = n, m = m, e = e, w = w, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
LS_CD_root / (LS_CD_root + get_a12(x = LS_CD_root, A = A, n = n)) # mutant growth rate = 0 as complete generalist

# Evolution favors incomplete specialization 
LS_CD_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = LS_CD_root, linetype = "dotted") +
  ylab("Mutant growth rate")

# We can also see this by looking at how small deviations from the incomplete specialist phenotype decrease mutant growth rate
C1_ESS_2C_2R_LawlorSmith(x=LS_CD_root, A = A, n = n, m = m, e = e, w = w, dx = 0.01)
C1_ESS_2C_2R_LawlorSmith(x=LS_CD_root, A = A, n = n, m = m, e = e, w = w, dx = -0.01)


## Linear tradeoff

n <- 1

# Plot
LS_LN_neg <- curve(C1_ESS_2C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, dx = -1e-10), 0 + 0.01, A - 0.01)
LS_LN_pos <- curve(C1_ESS_2C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, dx = 1e-10), 0 + 0.01, A - 0.01)
LS_LN_df <- data.frame(x = c(LS_LN_neg$x, LS_LN_pos$x), 
                        y = c(LS_LN_neg$y, LS_LN_pos$y), 
                        dx = c(rep(-1e-10, length(LS_LN_neg$x)), rep(1e-10, length(LS_LN_pos$x))))

# no root (at least for this value of w)
LS_LN_root <- uniroot(C1_ESS_2C_2R_LawlorSmith, A = A, n = n, m = m, e = e, w = w, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
LS_LN_root / (LS_LN_root + get_a12(x = LS_LN_root, A = A, n = n)) 

uniroot(C1_ESS_2C_2R_LawlorSmith, A = A, n = n, m = m, e = e, w = w, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))

# for this value of w, increasing a11 is always favored
LS_LN_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = LS_LN_root, linetype = "dotted") + # no root solution
  ylab("Mutant growth rate")

## since there is no root, we can't examine how small deviations from root affect mutant growth rate
#C1_ESS_2C_2R_LawlorSmith(x=LS_LN_root, A = A, n = n, m = m, e = e, dx = 0.01)
#C1_ESS_2C_2R_LawlorSmith(x=LS_LN_root, A = A, n = n, m = m, e = e, dx = -0.01)


## One consumer food webs

## Concave up tradeoff
n <- 0.85

## Resource abundances

LS_CU_oneC_R1 <- sapply(X = seq(0, 2, 0.1), FUN = function(x) eq.R1_1C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, r, K, R1, R2, C1))
LS_CU_oneC_R2 <- sapply(X = seq(0, 2, 0.1), FUN = function(x) eq.R2_1C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, r, K, R1, R2, C1))
LS_CU_oneC_R_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                            y = c(LS_CU_oneC_R1, LS_CU_oneC_R2), 
                            R = c(rep("R1", length(LS_CU_oneC_R1)), rep("R2", length(LS_CU_oneC_R2))))

## REWRITE ## # Indicative of complete specialization, since deviations from a complete generalist result in positive mutant growth rates
LS_CU_oneC_R_df %>%
  ggplot(., aes(x=x, y=y, color=R)) +
  geom_path() +
  ylab("Resource Abundance")

LS_CU_oneC_R_df %>%
  ggplot(., aes(x=(x/(x+get_a12(x=x, A=A, n=n))), y=y, color=R)) +
  geom_path() +
  ylab("Resource Abundance")

LS_CU_oneC_R_df %>%
  ggplot(., aes(x=(w*x/(w*x+(1-w)*get_a12(x=x, A=A, n=n))), y=y, color=R)) +
  geom_path() +
  ylab("Resource Abundance")

# Plot
LS_CU_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, r, K, R1, R2, C1, dx = -0.01))
LS_CU_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, r, K, R1, R2, C1, dx = 0.01))
LS_CU_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                             y = c(LS_CU_oneC_neg, LS_CU_oneC_pos), 
                             dx = c(rep(-0.01, length(LS_CU_oneC_neg)), rep(0.01, length(LS_CU_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
LS_CU_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = LS_CU_root, linetype = "dotted") + # root for two consumer DOES NOT predict root for one consumer
  ylab("Mutant growth rate")

LS_CU_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=(x/(x+get_a12(x=x, A=A, n=n))), y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = LS_CU_root, linetype = "dotted") + # root for two consumer DOES NOT predict root for one consumer
  ylab("Mutant growth rate")

LS_CU_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=(w*x/(w*x+(1-w)*get_a12(x=x, A=A, n=n))), y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = LS_CU_root, linetype = "dotted") + # root for two consumer DOES NOT predict root for one consumer
  ylab("Mutant growth rate")

## Concave down tradeoff
n <- 1.15

# Plot
LS_CD_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, r, K, R1, R2, C1, dx = -0.01))
LS_CD_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, r, K, R1, R2, C1, dx = 0.01))
LS_CD_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                             y = c(LS_CD_oneC_neg, LS_CD_oneC_pos), 
                             dx = c(rep(-0.01, length(LS_CD_oneC_neg)), rep(0.01, length(LS_CD_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
LS_CD_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = LS_CD_root, linetype = "dotted") + # root for two consumer DOES NOT predict root for one consumer
  ylab("Mutant growth rate")


## Linear tradeoff
n <- 1

# Plot
LS_LN_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, r, K, R1, R2, C1, dx = -0.01))
LS_LN_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_LawlorSmith(x, A = A, n = n, m = m, e = e, w = w, r, K, R1, R2, C1, dx = 0.01))
LS_LN_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                             y = c(LS_LN_oneC_neg, LS_LN_oneC_pos), 
                             dx = c(rep(-0.01, length(LS_LN_oneC_neg)), rep(0.01, length(LS_LN_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
LS_LN_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = A/2, linetype = "dotted") + # root for two consumer DOES NOT predict root for one consumer
  ylab("Mutant growth rate")

LS_LN_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=(x/(x+get_a12(x=x, A=A, n=n))), y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = McCann_CU_root / (McCann_CU_root + get_a12(x = McCann_CU_root, A = A, n = n)), linetype = "dotted") +
  ylab("Mutant growth rate")

LS_LN_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=(w*x/(w*x+(1-w)*get_a12(x=x, A=A, n=n))), y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = w*McCann_CU_root / (w*McCann_CU_root + (1-w)*get_a12(x = McCann_CU_root, A = A, n = n)), linetype = "dotted") +
  ylab("Mutant growth rate")


###################################################################################################################
#####  McCann 
###################################################################################################################

## Adding handling time
h <- 0.4

## Functions

C1_ESS_2C_2R_McCann <- function(x, A, n, m, e, w, h, dx){
  
  # calculate a12 in terms of a11 (x; tradeoff)
  a11 = x
  a12 = A*(1-(a11/A)^n)^(1/n)
  
  # calculate Resource abundances at equilibrium
  eq.R1 = 1/(w*a11 + (1-w)*a12) * m/(e - h*m)
  eq.R2 = 1/(w*a11 + (1-w)*a12) * m/(e - h*m)
  
  # calculate mutant attack rates
  a11m = x + dx
  a12m = A*(1-(a11m/A)^n)^(1/n)
  
  # calculate mutant growth rate of C1 when rare
  e*(((w * eq.R1)/(w * eq.R1 + (1-w) * eq.R2)) * a11m * eq.R1)/(1 + ((w * eq.R1)/(w * eq.R1 + (1-w) * eq.R2)) * a11m * h * eq.R1 + (((1-w) * eq.R2)/(w * eq.R1 + (1-w) * eq.R2)) * a12m * h * eq.R2) + e*((((1-w) * eq.R2)/(w * eq.R1 + (1-w) * eq.R2)) * a12m * eq.R2)/(1 + ((w * eq.R1)/(w * eq.R1 + (1-w) * eq.R2)) * a11m * h * eq.R1 + (((1-w) * eq.R2)/(w * eq.R1 + (1-w) * eq.R2)) * a12m * h * eq.R2) - m
}

C1_ESS_1C_2R_McCann <- function(x, A, n, m, e, w, h, dx, r, K, R1, R2, C1){
  
  # calculate a12 in terms of a11 (x; tradeoff)
  a11 <- x
  a12 <- A*(1-(a11/A)^n)^(1/n)
  
  # calculate Resource abundances at equilibrium
  get_steady <- safe.runsteady(y = c(R1 = R1, R2 = R2, C1 = C1), 
                               func = McCann_1C_2R, 
                               parms = c(a11 = a11, a12 = a12, 
                                         r1 = r, r2 = r,
                                         K1 = K, K2 = K,
                                         e11 = e, e12 = e,
                                         w11 = w, w12 = (1-w),
                                         h11 = h, h12 = h,
                                         m1 = m)) 
  
  eq.R1 = get_steady$y["R1"][[1]]
  eq.R2 = get_steady$y["R2"][[1]]
  
  # calculate mutant attack rates
  a11m = x + dx
  a12m = A*(1-(a11m/A)^n)^(1/n)
  
  # calculate mutant growth rate of C1 when rare
  e*(((w * eq.R1)/(w * eq.R1 + (1-w) * eq.R2)) * a11m * eq.R1)/(1 + ((w * eq.R1)/(w * eq.R1 + (1-w) * eq.R2)) * a11m * h * eq.R1 + (((1-w) * eq.R2)/(w * eq.R1 + (1-w) * eq.R2)) * a12m * h * eq.R2) + e*((((1-w) * eq.R2)/(w * eq.R1 + (1-w) * eq.R2)) * a12m * eq.R2)/(1 + ((w * eq.R1)/(w * eq.R1 + (1-w) * eq.R2)) * a11m * h * eq.R1 + (((1-w) * eq.R2)/(w * eq.R1 + (1-w) * eq.R2)) * a12m * h * eq.R2) - m
}

## Two consumer food webs

## Concave up tradeoff
n <- 0.85

# Plot
McCann_CU_neg <- curve(C1_ESS_2C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, dx = -1e-10), 0 + 0.01, A - 0.01)
McCann_CU_pos <- curve(C1_ESS_2C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, dx = 1e-10), 0 + 0.01, A - 0.01)
McCann_CU_df <- data.frame(x = c(McCann_CU_neg$x, McCann_CU_pos$x), 
                       y = c(McCann_CU_neg$y, McCann_CU_pos$y), 
                       dx = c(rep(-1e-10, length(McCann_CU_neg$x)), rep(1e-10, length(McCann_CU_pos$x))))

McCann_CU_root <- uniroot(C1_ESS_2C_2R_McCann, A = A, n = n, m = m, e = e, w = w, h = h, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
McCann_CU_root / (McCann_CU_root + get_a12(x = McCann_CU_root, A = A, n = n)) # mutant growth rate = 0 as complete generalist

## REWRITE ## # Indicative of complete specialization, since deviations from a complete generalist result in positive mutant growth rates
McCann_CU_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = McCann_CU_root, linetype = "dotted") +
  ylab("Mutant growth rate")

McCann_CU_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=(x/(x+get_a12(x=x, A=A, n=n))), y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = McCann_CU_root / (McCann_CU_root + get_a12(x = McCann_CU_root, A = A, n = n)), linetype = "dotted") +
  ylab("Mutant growth rate")

McCann_CU_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=(w*x/(w*x+(1-w)*get_a12(x=x, A=A, n=n))), y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = w*McCann_CU_root / (w*McCann_CU_root + (1-w)*get_a12(x = McCann_CU_root, A = A, n = n)), linetype = "dotted") +
  ylab("Mutant growth rate")

## REWRITE ## # We can also see this by looking at how small deviations from complete generalist phenotype increase mutant growth rate
C1_ESS_2C_2R_McCann(x=McCann_CU_root, A = A, n = n, m = m, e = e, w = w, h = h, dx = 0.01)
C1_ESS_2C_2R_McCann(x=McCann_CU_root, A = A, n = n, m = m, e = e, w = w, h = h, dx = -0.01)


## Concave down tradeoff

n <- 1.15

# Plot
McCann_CD_neg <- curve(C1_ESS_2C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, dx = -1e-10), 0 + 0.01, A - 0.01)
McCann_CD_pos <- curve(C1_ESS_2C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, dx = 1e-10), 0 + 0.01, A - 0.01)
McCann_CD_df <- data.frame(x = c(McCann_CD_neg$x, McCann_CD_pos$x), 
                       y = c(McCann_CD_neg$y, McCann_CD_pos$y), 
                       dx = c(rep(-1e-10, length(McCann_CD_neg$x)), rep(1e-10, length(McCann_CD_pos$x))))

McCann_CD_root <- uniroot(C1_ESS_2C_2R_McCann, A = A, n = n, m = m, e = e, w = w, h = h, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
McCann_CD_root / (McCann_CD_root + get_a12(x = McCann_CD_root, A = A, n = n)) # mutant growth rate = 0 as complete generalist

# Evolution favors incomplete specialization 
McCann_CD_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = McCann_CD_root, linetype = "dotted") +
  ylab("Mutant growth rate")

# We can also see this by looking at how small deviations from the incomplete specialist phenotype decrease mutant growth rate
C1_ESS_2C_2R_McCann(x=McCann_CD_root, A = A, n = n, m = m, e = e, w = w, h = h, dx = 0.01)
C1_ESS_2C_2R_McCann(x=McCann_CD_root, A = A, n = n, m = m, e = e, w = w, h = h, dx = -0.01)


## Linear tradeoff

n <- 1

# Plot
McCann_LN_neg <- curve(C1_ESS_2C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, dx = -1e-10), 0 + 0.01, A - 0.01)
McCann_LN_pos <- curve(C1_ESS_2C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, dx = 1e-10), 0 + 0.01, A - 0.01)
McCann_LN_df <- data.frame(x = c(McCann_LN_neg$x, McCann_LN_pos$x), 
                       y = c(McCann_LN_neg$y, McCann_LN_pos$y), 
                       dx = c(rep(-1e-10, length(McCann_LN_neg$x)), rep(1e-10, length(McCann_LN_pos$x))))

# no root (at least for this value of w)
McCann_LN_root <- uniroot(C1_ESS_2C_2R_McCann, A = A, n = n, m = m, e = e, w = w, h = h, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))$root
McCann_LN_root / (McCann_LN_root + get_a12(x = McCann_LN_root, A = A, n = n)) 

uniroot(C1_ESS_2C_2R_McCann, A = A, n = n, m = m, e = e, w = w, h = h, dx = 1e-10, interval = c(0 + 0.01, A - 0.01))

# for this value of w, increasing a11 is always favored
McCann_LN_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = McCann_LN_root, linetype = "dotted") + # no root solution
  ylab("Mutant growth rate")

## since there is no root, we can't examine how small deviations from root affect mutant growth rate
#C1_ESS_2C_2R_McCann(x=McCann_LN_root, A = A, n = n, m = m, e = e, dx = 0.01)
#C1_ESS_2C_2R_McCann(x=McCann_LN_root, A = A, n = n, m = m, e = e, dx = -0.01)


## One consumer food webs

## Concave up tradeoff
n <- 0.85

# Plot
McCann_CU_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, r, K, R1, R2, C1, dx = -0.01))
McCann_CU_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, r, K, R1, R2, C1, dx = 0.01))
McCann_CU_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                            y = c(McCann_CU_oneC_neg, McCann_CU_oneC_pos), 
                            dx = c(rep(-0.01, length(McCann_CU_oneC_neg)), rep(0.01, length(McCann_CU_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
McCann_CU_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = McCann_CU_root, linetype = "dotted") + # root for two consumer DOES NOT predict root for one consumer
  ylab("Mutant growth rate") 

## Concave down tradeoff
n <- 1.15

# Plot
McCann_CD_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, r, K, R1, R2, C1, dx = -0.01))
McCann_CD_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, r, K, R1, R2, C1, dx = 0.01))
McCann_CD_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                            y = c(McCann_CD_oneC_neg, McCann_CD_oneC_pos), 
                            dx = c(rep(-0.01, length(McCann_CD_oneC_neg)), rep(0.01, length(McCann_CD_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
McCann_CD_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = McCann_CD_root, linetype = "dotted") + # root for two consumer DOES NOT predict root for one consumer
  ylab("Mutant growth rate") 

McCann_CD_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=(x/(x+get_a12(x=x, A=A, n=n))), y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = McCann_CU_root / (McCann_CU_root + get_a12(x = McCann_CU_root, A = A, n = n)), linetype = "dotted") +
  ylab("Mutant growth rate")

McCann_CD_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=(w*x/(w*x+(1-w)*get_a12(x=x, A=A, n=n))), y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = w*McCann_CU_root / (w*McCann_CU_root + (1-w)*get_a12(x = McCann_CU_root, A = A, n = n)), linetype = "dotted") +
  ylab("Mutant growth rate")


## Linear tradeoff
n <- 1

# Plot
McCann_LN_oneC_neg <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, r, K, R1, R2, C1, dx = -0.01))
McCann_LN_oneC_pos <- sapply(X = seq(0, 2, 0.1), FUN = function(x) C1_ESS_1C_2R_McCann(x, A = A, n = n, m = m, e = e, w = w, h = h, r, K, R1, R2, C1, dx = 0.01))
McCann_LN_oneC_df <- data.frame(x = c(seq(0, 2, 0.1), seq(0, 2, 0.1)), 
                            y = c(McCann_LN_oneC_neg, McCann_LN_oneC_pos), 
                            dx = c(rep(-0.01, length(McCann_LN_oneC_neg)), rep(0.01, length(McCann_LN_oneC_pos))))

# Indicative of complete generalist, since deviations from a complete generalist result in negative mutant growth rates
McCann_LN_oneC_df %>%
  mutate(dx = as.factor(dx)) %>%
  ggplot(., aes(x=x, y=y, color=dx)) +
  geom_path() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #geom_vline(xintercept = A/2, linetype = "dotted") + # root for two consumer DOES NOT predict root for one consumer
  ylab("Mutant growth rate")+
  coord_cartesian(ylim = c(-0.01,0.01))
