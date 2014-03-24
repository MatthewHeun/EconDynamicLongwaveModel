require(deSolve)
require(FME)
require(lattice)


# A test function for R's DAE solvers.
# Example DAE system:
#   dy1/dt = a
#   y2 = y1 + b
#   y1_0 = 0
# The state variable is y1. 
# The algebraic variable is y2.
# The initial condition (need one for each state variable) is y1_0 = 0.
# In the required form as residuals:
#   dy1/dt - a = 0
#   y2 - y1 - b = 0
# And initial condition:
#   y1_0 = 0

# Let's define the parameters for the model first
parms <- c(a=5, b=2)

# Standard deviation for the data scatter
addNoise <- TRUE
sigma <- 1

# Specify the times at which we would like solutions reported
solveTimes <- seq(0, 2, 0.01) 
y1_init <- 10 # An initial condition to be used everywhere.

# Set up a function that solves the system using the DAE solver (daspk) in the deSolve package
# This function takes as arguments:
#   * times: a list of times at which you want the solution reported,
#   * parms: parameters to the function,
#   and either
#   * y1_init: an initial value for y1
#   * y_init: a vector of initial y values with names y1 and y2.
# Internally, this function calculates the initial value for the derivatives at time 0.
# I have buried the residuals function inside simpleSystemSolve, so that
# the code is nicely compartmentalized.
simpleDAE <- function(times, 
                      y1_init,
                      y_init=c(y1=y1_init, y2=as.vector(y1_init+parms["b"])),
                      parms){

  # Function that calculates residuals.
  # Bury inside simpleSystemSolve, because we don't need it anywhere else.
  residuals <- function(t, y, dydt, parms){
    # The block of code to calculate the residuals goes here.
    with(as.list(c(y, dydt, parms)), {
      # The "with" statement allows us to simplify the equations.
      # Because we're using "with," we can refer to the variables by their names
      # in the y, dydt, and parms numeric lists.
      # We don't need to de-reference them with statemens such as
      # a <- parms[["a"]]
      # Be careful for name collisions, though.
      R1 <- dy1dt - a
      R2 <- y2 - y1 - b
      # I'm adding dydt to the output so that we can see it later.
      out <- list(c(R1, R2))
      return(out)
    })
  }
  # Find a consistent set of derivatives and algebraic variables at the initial time.
  # For this problem, it is easy. By inspection, we can see what the derivatives are at the initial time.
  # But, for more complicated systems, we may need to execute a solve step.
  dydt_init <- c(dy1dt=as.vector(parms["a"]))
  # Use the residual function, the initial state, and the derivatives function
  # to integrate the DAE system forward in time with the daspk function.
  result <- daspk(y=y_init, 
                  times=solveTimes, 
                  dy=dydt_init,
                  parms=parms, 
                  res=residuals)
  out <- as.data.frame(result)
  return(out)
}

# Run the model using default values of the parameters (a and b)
unperturbedModel <- simpleDAE(times=solveTimes, y1_init=y1_init, parms=parms)

# If we know some "historical" data for y1 and y2, can we estimate values for a and b?

# Add some random noise to the solution to create the "historical" data that we will fit.
historical <- as.data.frame(solveTimes)
colnames(historical) <- c("time")
historical$y1 <- historical$time * parms["a"] + y1_init
historical$y2 <- historical$y1 + parms["b"] 
if (addNoise){
  historical$y1 <- historical$y1 + rnorm(sd=sigma, n=length(historical$time))
  historical$y2 <- historical$y2 + rnorm(sd=sigma, n=length(historical$time))
}

# Need to define a "cost" function.
# We'll use the modCost function in package FME to help with this.
# The model "cost" is a function of the parameters (parms)
simpleDAECost <- function(parms){
  pred <- simpleDAE(times=solveTimes, parms=parms, y1=y1_init)
  cost <- modCost(model=pred, obs=historical, x="time")
  return(cost)
}

# Try the cost function.
cost <- simpleDAECost(parms)
# We can plot the residuals as a function of time.
# We expect nice, random scatter for a good model.
par(mfrow=c(1, 1))
plot(cost, xlab="time")

# Now use DAECost in the modFit function to estimate values for the parameters a and b.
initGuess <- c(a=1.0, b=1.0)
model <- modFit(f=simpleDAECost, p=initGuess)
# The parameters that provide the best fit are in model$par.
# We can use those parameters in simpleDAE to find our fitted prediction.
fitted <- simpleDAE(times=solveTimes, y1_init=y1_init, parms=model$par)

# Make a plot of historical data and the fit 
par(mfrow=c(1, 2))
plot(y1~time, data=historical, xlab="time", ylab="y1")
lines(y1~time, data=fitted)
plot(y2~time, data=historical, xlab="time", ylab="y2")
lines(y2~time, data=fitted)
par(mfrow=c(1, 1))

