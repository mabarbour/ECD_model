library(tidyverse)
test <- read.csv("symmetric_displacement_data.csv")

select(test, Model, Competitor, Time, sequence, feasibility, n, a11:a22, w11:w22, R1:C2) %>%
  mutate(a11=round(a11,3),
         a12=round(a12,3),
         a21=round(a21,3),
         a22=round(a22,3),
         R1=round(R1,2),
         R2=round(R2,2),
         C1=round(C1,2),
         C2=round(C2,2)) %>%
  arrange(Model, n, Competitor, Time)

filter(test, Model=="MacArthur")
filter(test, Model=="LawlorSmith")
filter(test, Model=="McCann")
