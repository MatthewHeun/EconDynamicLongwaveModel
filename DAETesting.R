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

# Let's define the parameters first
parms <- c(a=5, b=2)

# Specify the times at which we would like solutions reported
solveTimes <- seq(0, 2, 0.1) 

# Set up a function that solves the system using the DAE solver (daspk) in the deSolve package
# This function takes a list of times at which you want the solution reported (times),
# a set of initial values (y_init), and parameters to the function (parms).
# Internally, it calculates the initial value for the derivatives at time 0.
# I have buried the residuals function inside simpleSystemSolve, so that
# the code is nicely compartmentalized.
simpleModel <- function(times, 
                        y_init=c(y1=0, y2=as.vector(y1+parms["b"])), 
                        parms=c(a=5, b=2)){

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
      res1 <- dy1dt - a
      res2 <- y2 - y1 - b
      out <- list(c(res1, res2), a=a, b=b)
      return(out)
    })
  }
  # Find a consistent set of derivatives at the initial time.
  # For this problem, it is easy. By inspection, we can what the derivatives are that the initial time.
  # But, for more complicated systems, we may need to execute a solve step.
  dydt_init <- c(dy1dt=as.vector(parms["a"]))
  # Use the residual function and the initial state and derivatives
  # to integrate the DAE system forward in time with the daspk function.
  result <- daspk(y=y_init, 
                  times=solveTimes, 
                  dy=dydt_init,
                  parms=parms, 
                  res=residuals)
  out <- as.data.frame(result)
  return(out)
}

# Run the model
result <- simpleModel(times=solveTimes)

# Show a graph of the results
print(xyplot(y1+y2~time, data=result))

# If we know some "historical" data for y1 and y2, can we estimate values for a and b?
# Add some random noise to the historical data.
historical <- data.frame(result$time)
colnames(historical) <- c("time")
historical$y1 <- seq(0, 10, 0.5) + rnorm(sd=0.1, n=length(historical$time))
historical$y2 <- seq(2, 12, 0.5) + rnorm(sd=0.1, n=length(historical$time))

# Need to define a "cost" function.
# We'll use the modCost function in package FME to help with this.
# The model "cost" is a function of the parameters (parms)
simpleModelCost <- function(parms){
  pred <- simpleModel(times=solveTimes, parms=parms)
  cost <- modCost(model=pred, obs=historical, x="time")
  return(cost)
}

# Try the cost function.
cost <- simpleModelCost(parms)
# The result will be different each time, due to randomness added to the historical data
print(cost$model) 

# Now use the cost function to estimate values for the parameters a and b.
fitted <- modFit(f=simpleModelCost, p=c(a=4.5, b=1.5))

