require(deSolve)
require(car)

# Read historical data
data <- read.table(file="data/USData.txt", header=TRUE, sep="\t")

# Direct calculations from historical data
data$L_Y <- data$L - data$L_A # workers
w_worker_0 <- data$w_worker[1] # Initial wage
Y_0  <- data$Y[1] # Initial GDP
N_0 <- data$N[1]   # Initial population
L_Y_0 <- data$L_Y[1] # Initial L_Y
data$y <- data$Y / Y_0 # indexed GDP
data$n <- data$N / N_0 # indexed population
data$l_Y <- data$L_Y / L_Y_0 # indexed labor
data$b <- data$B / data$N # Births/year-capita
data$d <- data$D / data$N # Deaths/year-capita
data$tau <- data$L / data$N # Fraction of workers in each year
data$c <- data$w_worker * data$tau # consumption in 2005$/year-capita. Same as wages in 2005$/year-capita.
data$nu <- data$b / (1.0 - data$tau)
data$pi <- data$L_A / data$L
data$wALA <- data$w_worker * data$L_A
data$piY <- data$pi * data$Y
data$zeta <- data$wALA / data$piY
b_bar = 0.0 # Asymptotic birth rate as consumption goes to infinity
data$b_tilde <- data$b - b_bar # Birth rate above asymptote, b_tilde
c_bar = 5655 # 2005$/year, Consumption for subsistence level
data$c_tilde <- data$c - c_bar # Consumption above subsistence level, c_tilde
d_bar <- 0 # Asymptotic mortality rate
data$z <- data$c / c_bar - 1

# Some initial values
y_0 <- data$y[1] # initial indexed GDP
n_0 <- data$n[1] # initial indexed population
b_0 <- data$b[1] # initial birth rate
d_0 <- data$d[1] # initial death rate
nu_0 <- data$nu[1] 
pi_0 <- data$pi[1]
zeta_0 <- data$zeta[1]
b_tilde_0 <- data$b_tilde[1]
c_tilde_0 <- data$c_tilde[1]
z_0 <- data$z[1]

# Estimate mortality parameters
mortalityModel <- d ~ 1.0 / (omega_1*z^omega_2 + (1.0 / (z_0*(d_0-d_bar)) - omega_1*z_0^(omega_2-1.0))*z) + d_bar
control <- nls.control(maxiter=200, 
                       tol=1e-05, 
                       minFactor=1/1024,
                       printEval=FALSE, #Tells whether to print details of curve fit process.
                       warnOnly=TRUE)
start <- list(omega_1=116 , omega_2=-0.63)
dModel <- nls(formula=mortalityModel, data=data, start=start, control=control)
omega_1 <- coef(dModel)["omega_1"]
omega_2 <- coef(dModel)["omega_2"]
omega_3 <- deltaMethod(dModel, "1/(z_0*(d_0-d_bar)) - omega_1*z_0^(omega_2-1.0)")[1, "Estimate"]
names(omega_3) <- "omega_3"

# Estimate utility model
utilityModel <- b_tilde ~ (w_worker/w_worker_0)^eta * (c_tilde/c_tilde_0)^(gamma/eta) * b_tilde_0
start <- list(eta=0.95, gamma=0.33)
# eta <- 0.95
# start <- list(gamma=0.33)
# uModel <- nls(formula=utilityModel, data=calculatedHistoricalData, start=start, control=control)
uModel <- nls(formula=utilityModel, data=data, start=start, control=control)
eta <- coef(uModel)["eta"]
gamma <- coef(uModel)["gamma"]
G <- w_worker_0 * b_tilde_0^eta / c_tilde_0^gamma


# 
# Steady state model
#
ssResid <- function(p, constraints, year){
  y <- p["y"]
  Y <- approx(x=data$Year, y=data$Y, xout=year)$y
  n <- p["n"]
  N <- approx(x=data$Year, y=data$N, xout=year)$y
  R1 <- as.vector(Y/Y_0 - y)
  R2 <- as.vector(N/N_0 - n)
  
  return(c(R1=R1, R2=R2))
}

ssYear <- 1985
ssParms <- c(dadt=0, dndt=0)
p_init <- c(y=0, n=0)
ssModel <- BBsolve(p=p_init, fn=ssResid, constraints=ssParms, year=ssYear)

#
# DAE Model
#
jonesDAE <- function(times, y_init, parms){
  
  # Function that calculates residuals.
  # Bury inside simpleSystemSolve, because we don't need it anywhere else.
  residuals <- function(t, y, dydt, parms){
    # The block of code to calculate the residuals goes here.
    with(as.list(c(y, dydt, parms)), {
      # The "with" statement allows us to simplify the equations.
      # Because we're using "with," we can refer to the variables by their names
      # in the y, dydt, and parms numeric lists.
      R1 <- dndt - (b - d + m)*n
      out <- list(c(R1))
      return(out)
    })
  }
  # Find a consistent set of derivatives and algebraic variables at the initial time.
  # For now, putting b and d in the parms vector. 
  # b and d will move to the y vector later.
  dndt_0 <- as.numeric((parms["b"] - parms["d"] + parms["m"])*y_init["n"])
  dydt_init <- c(dndt=dndt_0)

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

solveTimes <- seq(1980, 2011, 1)
y_init <- c(n=n_0) # Add b, d, and other variables when we expande beyond n.
parms <- c(m=0.0038, b=b_0, d=d_0) # m (migration), b (birth rate), and d (death rate) are placeholders
result <- jonesDAE(times=solveTimes, y_init=y_init, parms=parms)




