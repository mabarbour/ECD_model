
## UPDATE DETAILS OF SIMULATION RESULTS BELOW ##


## Description: Analyze effects of evolution and community context on character displacement and ecological stability

## source general parameters and functions ----
source('code/00_general_parameters.R')
source('code/00_numerical_fxns.R')

## load and manipulate required data ----
evol.sim.df.4 <- read.csv('data/evol.sim.df.4sp.csv') %>% 
  tbl_df() %>% select(-X) %>% 
  mutate(special.C1 = abs((a11/(a11 + a12)) - 0.5),
         special.C2 = abs((a22/(a21 + a22)) - 0.5),
         effective.C1 = a11*w11 + a12*(1-w11), # effective attack rate of consumer 1
         effective.C2 = a22*w22 + a21*(1-w22)) # effective attack rate of consumer 2

evol.sim.df.4.1evo <- read.csv('data/evol.sim.df.4sp.1evo.csv') %>% 
  tbl_df() %>% select(-X) %>% 
  mutate(special.C1 = abs((a11/(a11 + a12)) - 0.5),
         effective.C1 = a11*w11 + a12*(1-w11)) # effective attack rate of consumer 1

evol.sim.df.C1 <- read.csv('data/evol.sim.df.C1.3sp.csv') %>% 
  tbl_df() %>% select(-X) %>% 
  mutate(special.C1 = abs((a11/(a11 + a12)) - 0.5),
         effective.C1 = a11*w11 + a12*(1-w11)) # effective attack rate of consumer 1

evol.sim.df.C2 <- read.csv('data/evol.sim.df.C2.3sp.csv') %>% 
  tbl_df() %>% select(-X) %>% 
  mutate(special.C2 = abs((a22/(a21 + a22)) - 0.5),
         effective.C2 = a22*w22 + a21*(1-w22)) # effective attack rate of consumer 2

## Summarize simulation results ----

# 4 species - coevolution
summary.4 <- list()
for(i in 1:max(evol.sim.df.4$sim.number)){
  temp.df <- filter(evol.sim.df.4, sim.number == i)
  end.seq <- max(temp.df$sequence)
  
  summary.4[[i]] <- data.frame(
    sim.number = i,
    end.seq.4 = end.seq,
    R1.end.4 = as.numeric(temp.df[end.seq,"R1"]),
    R2.end.4 = as.numeric(temp.df[end.seq,"R2"]),
    effective.C1.end.4 = as.numeric(temp.df[end.seq,"effective.C1"]),
    effective.C2.end.4 = as.numeric(temp.df[end.seq,"effective.C2"]),
    special.C1.end.4 = as.numeric(temp.df[end.seq,"special.C1"]),
    special.C2.end.4 = as.numeric(temp.df[end.seq,"special.C2"]),
    feas.4sp.end = as.numeric(temp.df[end.seq,"feas.4sp"]),
    eigen.end.4 = as.numeric(temp.df[end.seq, "max.Re.eigen"]))
}

# 4 species - only C1 evolving
summary.4.C1evo <- list()
for(i in 1:max(evol.sim.df.4.1evo$sim.number)){
  temp.df <- filter(evol.sim.df.4.1evo, sim.number == i)
  end.seq <- max(temp.df$sequence)
  
  summary.4.C1evo[[i]] <- data.frame(
    sim.number = i,
    end.seq.4.C1evo = end.seq,
    R1.end.4.C1evo = as.numeric(temp.df[end.seq,"R1"]),
    R2.end.4.C1evo = as.numeric(temp.df[end.seq,"R2"]),
    effective.C1.end.4.C1evo = as.numeric(temp.df[end.seq,"effective.C1"]),
    special.C1.end.4.C1evo = as.numeric(temp.df[end.seq,"special.C1"]),
    feas.4sp.end.C1evo = as.numeric(temp.df[end.seq,"feas.4sp"]),
    eigen.end.4.C1evo = as.numeric(temp.df[end.seq, "max.Re.eigen"]))
}

# 3 species C1
summary.C1.3 <- list()
for(i in 1:max(evol.sim.df.C1$sim.number)){
  temp.df <- filter(evol.sim.df.C1, sim.number == i)
  end.seq <- max(temp.df$sequence)
  
  summary.C1.3[[i]] <- data.frame(
    sim.number = i,
    end.seq.C1.3 = end.seq,
    R1.end.C1.3 = as.numeric(temp.df[end.seq,"R1"]),
    R2.end.C1.3 = as.numeric(temp.df[end.seq,"R2"]),
    effective.C1.end.3 = as.numeric(temp.df[end.seq,"effective.C1"]),
    special.C1.end.3 = as.numeric(temp.df[end.seq,"special.C1"]),
    feas.C1.3sp.end = as.numeric(temp.df[end.seq,"feas.3sp"]),
    eigen.end.C1.3 = as.numeric(temp.df[end.seq,"C1.max.eigen"]))
}

# 3 species C2
summary.C2.3 <- list()
for(i in 1:max(evol.sim.df.C2$sim.number)){
  temp.df <- filter(evol.sim.df.C2, sim.number == i)
  end.seq <- max(temp.df$sequence)
  
  summary.C2.3[[i]] <- data.frame(
    sim.number = i,
    end.seq.C2.3 = end.seq,
    effective.C2.end.3 = as.numeric(temp.df[end.seq,"effective.C2"]),
    special.C2.end.3 = as.numeric(temp.df[end.seq,"special.C2"]),
    feas.C2.3sp.end = as.numeric(temp.df[end.seq,"feas.3sp"]),
    eigen.end.C2.3 = as.numeric(temp.df[end.seq,"C2.max.eigen"]))
}

## join data sets together and compute some more output ----
ECD.df <- ldply(summary.4) %>% 
  left_join(., ldply(summary.4.C1evo)) %>%
  left_join(., ldply(summary.C1.3)) %>% 
  left_join(., ldply(summary.C2.3)) %>% 
  mutate(ECD.C1 = ifelse((special.C1.end.4 - special.C1.end.3) > 0, "div", "con"), 
         ECD.C2 = ifelse((special.C2.end.4 - special.C2.end.3) > 0, "div", "con"), 
         coevo.C1 = ifelse((special.C1.end.4 - special.C1.end.4.C1evo) > 0, "div", "con"),
         comm.C1 = ifelse((special.C1.end.4.C1evo - special.C1.end.3) > 0, "div", "con"), 
         ecd.R1 = ifelse((R1.end.4 - R1.end.C1.3) > 0, "more", "less"),
         ecd.R2 = ifelse((R2.end.4 - R2.end.C1.3) > 0, "more", "less"),
         stab.4 = ifelse(eigen.end.4 > 0, "unstable", "stable"),
         stab.4.C1evo = ifelse(eigen.end.4.C1evo > 0, "unstable", "stable"),
         stab.C1.3 = ifelse(eigen.end.C1.3 > 0, "unstable", "stable"),
         eigen.ECD.C1 = ifelse((-1*eigen.end.4 - (-1*eigen.end.C1.3)) > 0, "stabilize", "destabilize"),
         eigen.ECD.C2 = ifelse((-1*eigen.end.4 - (-1*eigen.end.C2.3)) > 0, "stabilize", "destabilize"),
         eigen.coevo.C1 = ifelse((-1*eigen.end.4 - (-1*eigen.end.4.C1evo)) > 0, "stabilize", "destabilize"),
         eigen.comm.C1 = ifelse((-1*eigen.end.4.C1evo - (-1*eigen.end.C1.3)) > 0, "stabilize", "destabilize")) %>%
  tbl_df()

## explore results: only comparing C1 ----

## all simulations stopped short of maxium possible mutation steps ----
any(ECD.df$end.seq.4 == mut.reps) 
any(ECD.df$end.seq.4.C1evo == mut.reps) 
any(ECD.df$end.seq.C1.3 == mut.reps)

## How does community context and coevolution affect species exclusion? ----

# exclusion only happens in 4 species scenario (43% of the time)
with(ECD.df, table(feas.4sp.end, feas.C1.3sp.end))
with(ECD.df, prop.table(table(feas.4sp.end, feas.C1.3sp.end)))

# of this 43%, 22% is due to the addition of a fourth species
with(ECD.df, table(feas.4sp.end.C1evo, feas.C1.3sp.end, feas.4sp.end))
with(ECD.df, prop.table(table(feas.4sp.end.C1evo, feas.C1.3sp.end, feas.4sp.end)))

# of this 43%, 22% is due to coevolution
with(ECD.df, table(feas.4sp.end, feas.4sp.end.C1evo, feas.C1.3sp.end)) 
with(ECD.df, prop.table(table(feas.4sp.end, feas.4sp.end.C1evo, feas.C1.3sp.end)))

## How do community context and coevolution affect stability? ----
# restrict to feasible models
f.ECD.df <- ECD.df %>% filter(feas.4sp.end == 1, feas.4sp.end.C1evo == 1, feas.C1.3sp.end == 1) 

## Qualitative stability ----
# a 4 species system became unstable 23% of the time. 
# 20% of this was due to a 4 species system being qualitatively less stable than a 3 species system.
with(f.ECD.df, table(stab.4, stab.C1.3)) 
with(f.ECD.df, prop.table(table(stab.4, stab.C1.3))) 

# of this 20%, 9% became unstable due to the addition of a 4th species.
with(f.ECD.df, table(stab.4.C1evo, stab.C1.3, stab.4)) 
with(f.ECD.df, prop.table(table(stab.4.C1evo, stab.C1.3, stab.4)))

# of this 20%, 12% became unstable as a result of coevolution
with(f.ECD.df, table(stab.4, stab.4.C1evo, stab.C1.3)) 
with(f.ECD.df, prop.table(table(stab.4, stab.4.C1evo, stab.C1.3))) 

## Quantitative stability ----
# focus on feasible and stable solutions
fs.ECD.df <- ECD.df %>% filter(feas.4sp.end == 1, feas.4sp.end.C1evo == 1, feas.C1.3sp.end == 1, eigen.end.4 < 0, eigen.end.4.C1evo < 0, eigen.end.C1.3 < 0)

# 83% of simulations were less stable with 4 species and coevolution
with(fs.ECD.df, table(eigen.ECD.C1))
with(fs.ECD.df, prop.table(table(eigen.ECD.C1)))

# of this 83%, 68% was due to the addition of a 4th species. While 14% was due to coevolution.
with(fs.ECD.df, table(eigen.ECD.C1, eigen.comm.C1))
with(fs.ECD.df, prop.table(table(eigen.ECD.C1, eigen.comm.C1)))

# quantitatively, coevolution tends to destabilize rather than stabilize.
with(fs.ECD.df, table(eigen.coevo.C1))
with(fs.ECD.df, prop.table(table(eigen.coevo.C1)))

## How does community context and coevolution affect character displacement? ----

# consider re-running results by also looking at scenarios where C2 diverges as well
with(fs.ECD.df, table(ECD.C1, ECD.C2))
with(fs.ECD.df, prop.table(table(ECD.C1, ECD.C2)))

# divergent character displacement occurred 77% of the time (for C1) 
with(fs.ECD.df, table(ECD.C1))
with(fs.ECD.df, prop.table(table(ECD.C1)))

# of this 77%, the addition of a fourth species was responsible for 72% of the occassions. Coevolution only contributed to divergence 5% of the time. 
with(fs.ECD.df, table(ECD.C1, comm.C1))
with(fs.ECD.df, prop.table(table(ECD.C1, comm.C1)))

# Interestingly, coevolution resulted in convergence (relative to 4 sp. C1 evolving) 80%. THIS IS RATHER SURPRSING!
with(fs.ECD.df, table(coevo.C1)) 
with(fs.ECD.df, prop.table(table(coevo.C1)))

## How does character divergence affect community stability? ----

# focus on divergence
div.ECD.df <- filter(fs.ECD.df, ECD.C1 == "div", ECD.C2 == "div")

# of the simulations resulting in divergent character displacement. 85% of them became less stable than the 3 species system. 
with(div.ECD.df, table(eigen.ECD.C1, eigen.ECD.C2))
with(div.ECD.df, prop.table(table(eigen.ECD.C1, eigen.ECD.C2)))

# of this 85%, 66% can be attributed to the addition of a 4th species. Coevolution only contributed 19% of the time.
with(div.ECD.df, table(comm.C1, eigen.comm.C1, eigen.ECD.C1))
with(div.ECD.df, prop.table(table(comm.C1, eigen.comm.C1, eigen.ECD.C1)))

# interestingly, coevolution caused convergence/destabilization (relative to 4 sp. and C1 evolving) 49% of the time when there was divergent character displacement and net destabilization. 
with(div.ECD.df, table(coevo.C1, eigen.coevo.C1, eigen.ECD.C1))
with(div.ECD.df, prop.table(table(coevo.C1, eigen.coevo.C1, eigen.ECD.C1)))

## How does character divergence affect resource abundances? ----

# 73% of the time, divergent displacement results in lower abundances of both resources. 
with(div.ECD.df, table(ecd.R1, ecd.R2))
with(div.ECD.df, prop.table(table(ecd.R1, ecd.R2)))


#### below may be useful.... ----
# filter to only feasible solutions
feas.ECD.df <- ECD.df %>% filter(feas.4sp.end == 1, feas.C1.3sp.end == 1, feas.C2.3sp.end == 1)

# trait changes
with(feas.ECD.df, table(C1.3, C1.4)) # most common: specialize in 4, generalize in 3
with(feas.ECD.df, table(C2.3, C2.4)) # most common: specialize in 4, generalize in 3

with(feas.ECD.df, table(C1.4, C2.4)) # most common: both specializing
with(feas.ECD.df, table(C1.3, C2.3)) # most common: both generalizing

# all together, still see that the most common is for a consumer to specialize in 4, but generalize in 3.
with(feas.ECD.df, table(C1.3, C1.4, C2.3, C2.4))

## Destabilization is only slight more common in 4 species scenario
with(feas.ECD.df, table(eigen.4))
with(feas.ECD.df, table(eigen.C1.3))
with(feas.ECD.df, table(eigen.C2.3))

with(feas.ECD.df, hist(eigen.diff.C1.3))
plot(eigen.end.4 ~ eigen.end.C1.3, feas.ECD.df)

with(feas.ECD.df, table(eigen.4, eigen.C1.3))
with(feas.ECD.df, table(eigen.C1.3, eigen.C2.3))

## qualitative instability is about twice as common in 4 species scenario
with(feas.ECD.df, table(stab.4))
with(feas.ECD.df, table(stab.C1.3))
with(feas.ECD.df, table(stab.C2.3))

# filter some situations to see how common they are
feas.ECD.df %>% 
  filter(eigen.4 == "destabilize", eigen.C1.3 == "stabilize", eigen.C2.3 == "stabilize") %>% 
  select(C1.3, C1.4, C2.3, C2.4) %>%
  table()

feas.ECD.df %>% filter(C1.4 == "spec", C2.4 == "spec", C1.3 == "gen", C2.3 == "gen")  %>% 
with(comm.ECD.df, table(stab.4))
with(comm.ECD.df, table(eigen.4))

with(comm.ECD.df, table(stab.C1.3))
with(comm.ECD.df, table(eigen.C1.3))

uncomm.ECD.df <- feas.ECD.df %>% filter(C1.4 == "gen", C2.4 == "spec", C1.3 == "gen", C2.3 == "gen") 
with(uncomm.ECD.df, table(stab.4))
with(uncomm.ECD.df, table(eigen.4))

with(uncomm.ECD.df, table(stab.C1.3))
with(uncomm.ECD.df, table(eigen.C1.3))

##### below may be useful ----
# majority of feasible and stable solutions result in character divergence
ECD.df %>% 
  filter(feas.4sp.end == 1, feas.3sp.C1.end == 1, 
         stability.4 == "stable", stability.C1.3 == "stable") %>%
  select(ECD.C1, ECD.Stable.C1) %>%
  table()

ECD.df %>% 
  filter(feas.4sp.end == 1, feas.3sp.C2.end == 1, 
         stability.4 == "stable", stability.C2.3 == "stable") %>%
  select(ECD.C2, ECD.Stable.C2) %>%
  table()

# stability always decreased in 4 species communities and resulted in an unstable system about 23% of the time
# stability tended to decrease in 3 species communities, but increased in stability about 25% of the time
ECD.df %>% 
  filter(feas.4sp.end == 1, feas.3sp.C1.end == 1) %>%
  select(stability.4, stability.C1.3) %>%
  table() # so there are some situations where the 3 species system can become unstable whereas a 4 species system remains stable...

ECD.df %>% 
  filter(feas.4sp.end == 1, feas.3sp.C2.end == 1) %>%
  select(stability.4, stability.C2.3) %>%
  table() # same qualitative pattern

ECD.df %>% 
  filter(feas.4sp.end == 1, feas.3sp.C1.end == 1) %>%
  select(eigen.diff.4, eigen.diff.C1.3) %>%
  table()

ECD.df %>% 
  filter(feas.4sp.end == 1, feas.3sp.C1.end == 1) %>%
  select(eigen.diff.4, eigen.diff.C2.3) %>%
  table()

### below is old, ideas may be useful... ----
## Do we observe ecological character displacement?
hist(filter(ECD.df, steady.4sp.end == 1, steady.C1.3sp.end == 1)$C1.diff, breaks = 20)
hist(filter(ECD.df, steady.4sp.end == 1, steady.C2.3sp.end == 1)$C2.diff, breaks = 20)

## Does community context qualitative affect how evolution influences ecological stability?
table(with(ECD.df, steady.4sp.end - steady.C1.3sp.end), useNA = "always")
15/(15+23) # in 40% of the trials, evolution in 4 species systems decreases stability

table(with(ECD.df, steady.4sp.end - steady.C2.3sp.end), useNA = "always")
13/(13+25+1) # in 33% of the trials, evolution in 4 species systems decreases stability

## In the 4 species system, if it remained stable, did evolution decrease stability
hist(filter(ECD.df, steady.4sp.end == 1)$diff.eigen, breaks = 20)
hist(ECD.df$diff.eigen, breaks = 20)

plot(I(C1.4sp.end + C2.4sp.end) ~ steady.4sp.end, ECD.df) # tendency for qualitative stability when specialization is high...

table(ECD.df$steady.4sp.end, useNA = "always")
table(ECD.df$steady.C1.3sp.end, useNA = "always")
table(ECD.df$steady.C2.3sp.end, useNA = "always")

plot(C2.4sp.end ~ steady.4sp.end, ECD.df)

## Explore representative simulations ----
head(ECD.df)

## choose dataset
explore.df <- filter(evol.sim.df.C1, sim.number == 2)
explore.df.mat <- as.matrix(explore.df)
seq.num <- 27

## look at the important data
explore.df.mat[1:30,c("a11","a12","R1","R2","C1","C1.max.eigen","steady.3sp")]

## Check to make sure values matchup
runsteady(y = explore.df.mat[seq.num, 13:15], func = ECD_model.C1.3sp, parms = explore.df.mat[seq.num, 1:12], stol = 1e-4)
matplot.deSolve(ode(y = explore.df.mat[seq.num, 13:15], func = ECD_model.C1.3sp, parms = explore.df.mat[seq.num, 1:12], times = 1:1000))
eq.jac <- jacobian.full(y = explore.df.mat[seq.num, 13:15], func = ECD_model.C1.3sp, parms = explore.df.mat[seq.num, 1:12])
max(Re(eigen(eq.jac)$values))

##################### >Maybe useful below...

## only feasible solutions examining
#filter(ECD.df, stability.4sp == "exclusion", stability.C2.3sp == "stable", stability.C1.3sp == "stable") # no clear explanation yet. Maybe there is substantial difference in resource ratio?

filter(ECD.df, stability.4sp == "stable", stability.C1.3sp == "stable", stability.C2.3sp == "stable") # net divergence #stability.C2.3sp == "stable", 
filter(ECD.df, stability.4sp == "unstable", stability.C1.3sp == "stable") # net divergence #stability.C2.3sp == "stable", 

filter(ECD.df, stability.4sp == "exclusion", stability.C2.3sp == "exclusion", stability.C1.3sp == "exclusion") # interesting, at least 1 species is converging more in 4 species system.
filter(ECD.df, stability.4sp == "exclusion", stability.C2.3sp == "exclusion", stability.C1.3sp == "stable") # interesting, stable species always diverges in 4 species system.
filter(ECD.df, stability.4sp == "exclusion", stability.C2.3sp == "stable", stability.C1.3sp == "exclusion") # hmmmm, stable species isn't always diverging in 4 species system.

# interesting, evolution in a 4 species system can promote coexistence compared to evolution in 3 species system. Also, divergent ECD has occurred in the 4 species system in all of these cases.
filter(ECD.df, stability.4sp == "stable", stability.C2.3sp == "exclusion", stability.C1.3sp == "exclusion") 
filter(ECD.df, stability.4sp == "stable", stability.C2.3sp == "stable", stability.C1.3sp == "exclusion")
filter(ECD.df, stability.4sp == "stable", stability.C2.3sp == "exclusion", stability.C1.3sp == "stable")
filter(ECD.df, stability.4sp == "stable", stability.C2.3sp == "stable", stability.C1.3sp == "stable") # net divergence


## POSITIVE STABILITY, MIGHT ONLY BE WHEN AT LEAST ONE OF THE CONSUMERS HAS REALLY LOW INTERACTION STRENGTH (I.E. IN NON-EXCITABLE REGION OF PARAMETER SPACE!!!!)

## POTENTIALLY IMPORTANT QUESTION: WHICH REGIONS OF PARAMETER SPACE RESULT IN TRAJECTORIES TOWARD INSTABILITY VS. STABILITY?
# DO ALL OF THE ONES THAT TREND TOWARD INSTABILITY ULTIMATELY RESULT IN AN UNSTABLE SYSTEM?

output <- list()
for(i in 1:max(evol.sim.df.4$sim.number)){
  
  temp.df <- filter(evol.sim.df.4, sim.number == i)
  
  steps.to.exclude <- which(temp.df$steady == 0.5)[1]
  steps.to.cycle <- which(temp.df$max.eigen > 0)[1]
  
  # determine whether only one or both species successfully mutated
  mut.opp <- ifelse(length(with(filter(temp.df, mut.suc == 1), table(sp))) == 1, "one", "multiple") #ifelse(length(table(temp.df$sp)) == 1, "one", "multiple")
  
  # mutation constraint
  mut.const <- ifelse(any(temp.df$p.try == -1, na.rm = TRUE) == TRUE, "constraint", "none")
  
  # calculate increase/decrease in effective attack rates of both species
  diff.eff.C1 <- temp.df$effective.C1[dim(temp.df)[1]] - temp.df$effective.C1[1]
  diff.eff.C2 <- temp.df$effective.C2[dim(temp.df)[1]] - temp.df$effective.C2[1]
  max.eff.C1 <- temp.df$effective.C1[dim(temp.df)[1]]
  max.eff.C2 <- temp.df$effective.C2[dim(temp.df)[1]]
  
  # calculate increase/decrease in specialization of both species. 
  # Note that choice of following R1 or R2 doesn't matter since special.C1R1 + special.C1R2 = 1. In other words, we know special.C1R2 by knowing special.C1R1
  start.C1 <- abs(temp.df$special.C1R1[1] - 0.5)
  end.C1 <- abs(temp.df$special.C1R1[dim(temp.df)[1]] - 0.5)
  diff.C1 <- end.C1 - start.C1
  
  start.C2 <- abs(temp.df$special.C2R1[1] - 0.5)
  end.C2 <- abs(temp.df$special.C2R1[dim(temp.df)[1]] - 0.5)
  diff.C2 <- end.C2 - start.C2
  
  net.diff <- diff.C1 + diff.C2
  net.diff.qual <- ifelse(net.diff > 0, "net divergence", "net convergence")
  
  # calculate type of character displacement
  displacement.type <- ifelse(diff.C1 > 0 & diff.C2 > 0, "diverge",
                              ifelse(diff.C1 < 0 & diff.C2 < 0, "converge",
                                     ifelse(diff.C1 > 0 & diff.C2 < 0, "parallel",
                                            ifelse(diff.C1 < 0 & diff.C2 > 0, "parallel",
                                                   "other"))))
  

  # determine whether evolution results in unstable system or competitive exclusion
  # if qualitatively stable, calculate stability based on eigenvalue
  #if(any(temp.df$max.eigen > 0) == TRUE){
   # stability <- "unstable"
  #} 
  if(any(temp.df$steady == 0.5) == TRUE){
    stability <- "exclusion"
  }
  if(any(temp.df$steady < 1) == FALSE){
    stability <- "stable"
  }
  
  stable.region <- filter(temp.df, steady == 1)
  
  diff.stability <- -1*temp.df[dim(stable.region)[1], "max.eigen"] - # end.stability -
    -1*stable.region[1, "max.eigen"] # start.stability
  names(diff.stability) <- "diff.stability"
  
  stability.trend <- as.character(ifelse(diff.stability < 0, "destabilizing", ifelse(diff.stability > 0, "stabilizing", "NA")))
  
  # output
  output[[i]] <- data.frame(sim.number = i, steps.to.exclude = steps.to.exclude, steps.to.cycle = steps.to.cycle, diff.eff.C1 = diff.eff.C1, diff.eff.C2 = diff.eff.C2, max.eff.C1 = max.eff.C1, max.eff.C2 = max.eff.C2, diff.C1, diff.C2, net.diff, net.diff.qual, displacement.type = displacement.type, mut.opp = mut.opp, mut.const = mut.const, stability = stability, diff.stability = diff.stability, stability.trend = stability.trend)
}
output.df <- ldply(output)

filter(output.df, stability == "stable", stability.trend == "stabilizing")
filter(output.df, stability == "stable", stability.trend == "destabilizing", steps.to.cycle > 0)

#output.diff <- list()
#for(i in 1:max(evol.sim.df.4$sim.number)){
  
 # temp.df <- filter(evol.sim.df.4, sim.number == i)
  
  #output.diff.temp <- list()
  #for(j in 2:dim(temp.df)[1]){
   # diff.C1 <- abs(temp.df$special.C1R1[j-1] - 0.5) - abs(temp.df$special.C1R1[j] - 0.5)
    #diff.C2 <- abs(temp.df$special.C2R1[j-1] - 0.5) - abs(temp.df$special.C2R1[j] - 0.5)
    
    #diff.stability <- -1*temp.df$max.eigen[j-1] - (-1*temp.df$max.eigen[j]) 
    
    #diff.C <- max(c(diff.C1, diff.C2))
    
    #output.diff.temp[[j]] <- c(diff.C = diff.C, diff.stability = diff.stability)
  #}
  #output.diff[[i]] <- ldply(output.diff.temp)
#}
#output.diff.df <- ldply(output.diff)
#head(output.diff.df)
#plot(diff.stability ~ diff.C, filter(output.diff.df, diff.C > 0))

filter(output.df, net.diff.qual == "net convergence")

with(filter(output.df, mut.opp == "multiple"), table(displacement.type, stability.trend, mut.const))

## Does competitor coevolution result in divergence, convergence, or parallel displacement in foraging traits?
table(output.df$displacement.type)
table(output.df$net.diff.qual)

with(output.df, table(net.diff.qual, mut.opp))
with(output.df, table(net.diff.qual, mut.const))
with(output.df, table(mut.opp, stability))

## Does competitor coevolution result in qualitative changes in stability?
table(output.df$stability)

## Are different displacement types associated with differences in stability?
with(output.df, table(stability, displacement.type))
with(output.df, table(stability, net.diff.qual, mut.const))

with(output.df, table(stability.trend, displacement.type))

## Is the "other" displacement type always associated with no mutational opportunity?
with(output.df, table(displacement.type, mut.opp))

## For stable systems, does coevolution quantitatively affect stability?
table(output.df$stability.trend)

which(output.df$stability == "exclusion")
which(output.df$stability.trend == "stabilizing")
which(output.df$stability.trend == "destabilizing" & output.df$displacement.type == "parallel")

as.matrix(filter(evol.sim.df.4, sim.number == 1))[1:10, ] # positive stability scenario occurs with a consumer that has very low attack rates to begin with. Need to confirm whether this scenario was because the eigenvalue was not complex at that time.

plot(diff.stability ~ net.diff, output.df)
output.df$net.diff
output.df$diff.stability
which(output.df$diff.stability > 0)

## Explore representative simulations ----
which(output.df$stability == "stable" & output.df$displacement.type == "diverge")

## choose dataset
explore.df <- filter(evol.sim.df.4, sim.number == 2)
explore.df.mat <- as.matrix(explore.df)
seq.num <- dim(explore.df)[1]

## look at the important data
explore.df.mat[1:50,c("a11","a12","a21","a22","R1","R2","C1","C2","max.Re.eigen","steady.4sp")]

## Check to make sure values matchup
runsteady(y = explore.df.mat[seq.num, 21:24], func = ECD_model, parms = explore.df.mat[seq.num, 1:20])
stode(y = explore.df.mat[seq.num, 21:24], func = ECD_model, parms = explore.df.mat[seq.num, 1:20])
ode(y = explore.df.mat[seq.num, 21:24], func = ECD_model, parms = explore.df.mat[seq.num, 1:20], times = 1:1000)
eq.jac <- jac.norm.calc(data = c(explore.df.mat[seq.num, 1:20], explore.df.mat[seq.num, 21:24]), jac.mat = jac)
max(Re(eigen(eq.jac[1:4,1:4])$values))

## Dynamics matplot to check things
dynamic_matplot(param.vector = explore.df.mat[seq.num, 1:20],
                init.state = explore.df.mat[seq.num, 21:24],#c(R1 = 1,R2 = 1,C1 = 1, C2 = 1),#
                sim.length = 1000,
                model = ECD_model, 
                ylim = c(0,5),
                col = c("steelblue","steelblue","black","black"),
                lty = c(1,2,1,2),
                ylab = "Population Density", xlab = "Time")

## stability plot for one simulation
plot(-1*max.Re.eigen ~ sequence, explore.df)

## Create a plot of how mutations proceed over time 
# Using shaded arrows to draw the eye into the direction of mutations offer time.
# consider turning into ggplot
plot(0:1, 0:1, type = "n", 
     xlab = "Attack Rate \nSpecialization on R1", 
     ylab = "Attack Rate \nSpecialization on R2")
arrows(explore.df$special.C1R1[explore.df$sequence], explore.df$special.C1R2[explore.df$sequence],
       explore.df$special.C1R1[explore.df$sequence+1], explore.df$special.C1R2[explore.df$sequence+1], col = gray.colors(seq.num, start = 0.9, end = 0.1), length = 0.1, lwd = 5)
arrows(explore.df$special.C2R1[explore.df$sequence], explore.df$special.C2R2[explore.df$sequence],
       explore.df$special.C2R1[explore.df$sequence+1], explore.df$special.C2R2[explore.df$sequence+1], col = gray.colors(seq.num, start = 0.9, end = 0.1), length = 0.1, lwd = 5)
points(x = explore.df$special.C1R1[1], y = explore.df$special.C1R2[1], col = "black", bg = "white", pch = 21, cex = 5)
points(x = explore.df$special.C1R1[seq.num], y = explore.df$special.C1R2[seq.num], col = "black", bg = "white", pch = 21, cex = 5)
points(x = explore.df$special.C2R1[1], y = explore.df$special.C2R2[1], col = "black", bg = "white", pch = 21, cex = 5)
points(x = explore.df$special.C2R1[seq.num], y = explore.df$special.C2R2[seq.num], col = "black", bg = "white", pch = 21, cex = 5)
text(x = explore.df$special.C1R1[1], y = explore.df$special.C1R2[1], labels = "C1")
text(x = explore.df$special.C1R1[seq.num], y = explore.df$special.C1R2[seq.num], labels = "C1*")
text(x = explore.df$special.C2R1[1], y = explore.df$special.C2R2[1], labels = "C2")
text(x = explore.df$special.C2R1[seq.num], y = explore.df$special.C2R2[seq.num], labels = "C2*")

