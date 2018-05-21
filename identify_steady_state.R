## Function to identify steady-state equilibriums

# all runsteady simulations go for 1000 time steps with a steady state tolerance of 1e-4. 
runsteady.ps <- function(y, func = func, parms = parms, times = c(0,1000), stol = 1e-4){
  runsteady(y = y, times = times, func = func, parms = parms, stol = stol)
}

# using failwith() to make function give output even with a fatal error
safe.runsteady <- failwith(default = structure(list(y = c(R1 = NA, R2 = NA, C1 = NA, C2 = NA)), steady = 0), runsteady.ps)
