require(deSolve)

# A test function for R's DAE solvers.
# Example DAE system:
#   dy1/dt = a
#   y2 = y1 + b
#   y1_0 = 0
# In the required form:
#   dy1/dt - a = 0
#   y2 - y1 - b = 0
#   y1_0 = 0



y_init <- c(y1=0, y2=2) # Initial state 
dydt_init <- c(y1=5)    # Initial derivatives
solveTimes <- seq(0, 2, 0.1) # Times at which output is desired.
params <- c(a=5, b=2)

# Function that calculates residuals
residuals <- function(t, y, dy, parms){
  a <- parms[1]
  b <- parms[2]
  dy1dt <- dy[1]
  y1 <- y[1]
  y2 <- y[2]
  res1 <- dy1dt - a
  res2 <- y2 - y1 - b
  out <- list(c(res1, res2), a, b)
  return(out)
}

result <- daspk(y=y_init, times=solveTimes, dy=dydt_init,
                parms=params, res=residuals)
df <- data.frame(result)