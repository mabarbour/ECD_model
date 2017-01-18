
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

#evol.sim.df.4.1evo <- read.csv('data/evol.sim.df.4sp.1evo.csv') %>% 
#  tbl_df() %>% select(-X) %>% 
#  mutate(special.C1 = abs((a11/(a11 + a12)) - 0.5),
#         effective.C1 = a11*w11 + a12*(1-w11)) # effective attack rate of consumer 1

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
#summary.4.C1evo <- list()
#for(i in 1:max(evol.sim.df.4.1evo$sim.number)){
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
    R1.end.C2.3 = as.numeric(temp.df[end.seq,"R1"]),
    R2.end.C2.3 = as.numeric(temp.df[end.seq,"R2"]),
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
         ecd.C1.Rs = ifelse((R1.end.4 + R2.end.4 - R1.end.C1.3 - R2.end.C1.3) > 0, "more", "less"),
         ecd.C2.Rs = ifelse((R1.end.4 + R2.end.4 - R1.end.C2.3 - R2.end.C2.3) > 0, "more", "less"),
         stab.4 = ifelse(eigen.end.4 > 0, "unstable", "stable"),
         stab.C1.3 = ifelse(eigen.end.C1.3 > 0, "unstable", "stable"),
         stab.C2.3 = ifelse(eigen.end.C2.3 > 0, "unstable", "stable"),
         eigen.ECD.C1 = ifelse((-1*eigen.end.4 - (-1*eigen.end.C1.3)) > 0, "stabilize", "destabilize"),
         eigen.ECD.C2 = ifelse((-1*eigen.end.4 - (-1*eigen.end.C2.3)) > 0, "stabilize", "destabilize")) %>% 
  tbl_df()

## all simulations stopped short of maxium possible mutation steps
# confirms that simulations either reached a steady state, went extinct, or became unstable
any(ECD.df$end.seq.4 == mut.reps) 
any(ECD.df$end.seq.4.C1evo == mut.reps) 
any(ECD.df$end.seq.C1.3 == mut.reps)

## feasible data set
f.ECD.df <- ECD.df %>% 
  filter(feas.4sp.end == 1, feas.C1.3sp.end == 1, feas.C2.3sp.end == 1)

## feasible and stable data set
# most relevant for determining character displacement, which assumes an ESS has been reached
fs.ECD.df <- f.ECD.df %>% 
  filter(eigen.end.4 < 0, eigen.end.C1.3 < 0, eigen.end.C2.3 < 0) 

## How often did we observe character displacement? ----

# divergent character displacement occured the majority of time (53%)
with(fs.ECD.df, table(ECD.C1, ECD.C2))
with(fs.ECD.df, prop.table(table(ECD.C1, ECD.C2)))

## How does character divergence affect community stability? ----

# focus on divergent and stable solutions
div.ECD.df <- fs.ECD.df %>% filter(ECD.C1 == "div", ECD.C2 == "div")

# of the simulations resulting in divergent character displacement. 72% of them became less stable than the 3 species system. 
with(div.ECD.df, table(eigen.ECD.C1, eigen.ECD.C2))
with(div.ECD.df, prop.table(table(eigen.ECD.C1, eigen.ECD.C2)))

# focus on divergence and feasible (not necessarily stable) solutions
div.f.ECD.df <- f.ECD.df %>% filter(ECD.C1 == "div", ECD.C2 == "div")

# 4 species simulation resulted in qualitative instability 20% of time compared to stable 3 species systems
with(div.f.ECD.df, table(stab.C1.3, stab.C2.3, stab.4)) 
with(div.f.ECD.df, prop.table(table(stab.C1.3, stab.C2.3, stab.4))) 

## How does character divergence affect resource abundances? ----

# 83% of the time, divergent displacement results in lower total abundances of resources. 
with(div.ECD.df, table(ecd.C1.Rs, ecd.C2.Rs))
with(div.ECD.df, prop.table(table(ecd.C1.Rs, ecd.C2.Rs)))

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

