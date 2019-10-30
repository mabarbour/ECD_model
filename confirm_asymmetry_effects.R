test <- read.csv("asymmetric_displacement_data.csv")

# inspecting the data shows that the starting values are being appropriately switched. Also, we see the same qualitative pattern of displacement in the species that starts as a generalist (becomes a specialist).
select(test, Model, Competitor, sequence, Time, n, a11:a22, w11:w22, R1:C2, max.Re.eigen, max.Im.eigen) %>%
  mutate(a11=round(a11,3),
         a12=round(a12,3),
         a21=round(a21,3),
         a22=round(a22,3)) %>%
  arrange(Model, n, Competitor, Time)
