require(deSolve)
require(BB)
require(car)

# Read historical data
data <- read.table(file="data/USData.txt", header=TRUE, sep="\t")

# Direct calculations from historical data
data$L_Y <- data$L - data$L_A # workers
w_worker_0 <- data$w_worker[1] # Initial wage
Y_0  <- data$Y[1] # Initial GDP
N_0 <- data$N[1]   # Initial population
L_0 <- data$L[1]   # Initial labor force
L_Y_0 <- data$L_Y[1] # Initial L_Y
L_A_0 <- data$L_A[1] # Initial L_A
data$y <- data$Y / Y_0 # indexed GDP
data$n <- data$N / N_0 # indexed population
data$l_Y <- data$L_Y / L_Y_0 # indexed labor
data$b <- data$B / data$N # Births/year-capita
data$d <- data$D / data$N # Deaths/year-capita
data$tau <- data$L / data$N # Fraction of workers in each year
data$c <- data$w_worker * data$tau # consumption in 2005$/year-capita. Same as wages in 2005$/year-capita.
data$nu <- data$b / (1.0 - data$tau)
data$pi <- data$L_A / data$L
data$w_A <- data$w_worker # Jones says: "These wages will be equated in equilibrium by the free flow of labor between the two sectors."
data$wALA <- data$w_A * data$L_A # Question: should this be w_A * L_A? Do we have w_A?
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
l_Y_0 <- data$l_Y[1] # initial labor force rate
l_0 <- data$l[1]
nu_0 <- data$nu[1] 
pi_0 <- data$pi[1]
zeta_0 <- data$zeta[1]
b_tilde_0 <- data$b_tilde[1]
c_tilde_0 <- data$c_tilde[1]
z_0 <- data$z[1]

# Land model
t <- 1 # indexed land (never changes)

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
uModel <- nls(formula=utilityModel, data=data, start=start, control=control)
eta <- coef(uModel)["eta"]
gamma <- coef(uModel)["gamma"]
# G_0 <- w_worker_0 * b_tilde_0^eta / c_tilde_0^gamma


#' Residuals for the Jones model
#' 
#' Returns values for residuals for the Jones model. Note that \code{p} 
#' and \code{constraints} should NOT contain any repeated variable names.
#' 
#' @param p a vector with named values at which residuals are to be evaluated.
#' Both names and values in \coe{p} are important. 
#' The names in \code{p} are used internally to make the residual equations more readable.
#' The names in \code{p} are the variables to be solved for.
#' The residuals will be evaluated with the values in \code{p}.
#' @param constraints a vector of additional named values that constrain the solution.
ssResid <- function(p, constraints){
  # Ensure that p and constraints have no names in common.
  # If any names are in common, that is an error! You can't both solve for a variable (by putting it in the p vector) 
  # and use it as a constraint (by putting it in the constraints vector).
  if (length(commonNames <- intersect(names(p), names(constraints))) > 0){
    stop(paste("p and constraints both have variables named: "), commonNames)
  }
  with(as.list(c(p, constraints)), {
    Y <- approx(x=data$Year, y=data$Y, xout=year)$y # Interpolate to find Y at year
    N <- approx(x=data$Year, y=data$N, xout=year)$y # Interpolate to find N at year
    L_Y <- approx(x=data$Year, y=data$L_Y, xout=year)$y #Interpolate to find L_Y at year
    L <- approx(x=data$Year, y=data$L, xout=year)$y #Interpolate to find L at year
    w_worker <- approx(x=data$Year, y=data$w_worker, xout=year)$y #Interpolate to find w_worker at year
    w_A <- approx(x=data$Year, y=data$w_A, xout=year)$y #Interpolate to find w_A at year
    zeta <- approx(x=data$Year, y=data$zeta, xout=year)$y #Interpolate to find zeta at year
    nu <- approx(x=data$Year, y=data$nu, xout=year)$y #Interpolate to find nu at year

    R1 <- as.vector(Y/Y_0 - y) # y = Y/Y_0   as.vector() strips off the name
    R2 <- as.vector(N/N_0 - n) # n = N/N_0
    R3 <- as.vector(L_Y/L_Y_0 - l_Y) # l_Y = L_Y/L_Y_0
    R4 <- as.vector(L/L_0 - l) #l = L/L_0
    R5 <- as.vector(L/N - tau) #tau = L/N
    R6 <- as.vector(L - L_Y - L_A) #L_A = L - L_Y
    R7 <- as.vector(L_A / L_A_0 - l_A) #l_A = L_A/L_A_0
    R8 <- as.vector(w_worker * tau - c) #c = w_worker*tau
    R9 <- as.vector(w_A*L_A/(Y*zeta) - pi) #pi = w_A*L_A / (zeta*Y)
    R10 <- as.vector(c - c_bar - c_tilde) #c_tilde = c - c_bar
    R11 <- as.vector(nu*(1-tau) - b) #b = nu*(1-tau)
    R12 <- as.vector(b - b_bar - b_tilde) #b_tilde = b - b_bar
    R13 <- as.vector(c/c_bar - 1.0 - z) #z = c/c_bar - 1
    # We were working on this residual. We know everything in this equation except for d.
    # So, this equation should solve for d. But, we weren't getting a good solution.
#     R14 <- as.vector(1.0/(omega_1*z^omega_2 + omega_3*z) + d_bar - d) #d = 1/(omega_1*z^omega_2 + omega_3*z) + d_bar
    return(c(R1=R1, R2=R2, R3=R3, R4=R4, R5=R5, R6=R6, R7=R7, R8=R8, R9=R9, R10=R10, 
             R11=R11, R12=R12, R13=R13))
  })
}

p_init <- c(y=0, n=0, l_Y=0, l=0, tau=0, L_A=0, l_A=0, 
            c=0, pi=0, c_tilde=0, b=0, b_tilde=0, z=2) # Initial guess for the parameters that will be solved
ssParms <- c(dadt=0, dndt=0, year=1980, 
             omega_1=omega_1, omega_2=omega_2, omega_3=omega_3,
             d_bar=d_bar) # Constraint parameters for the model
ssModel <- BBsolve(p=p_init, fn=ssResid, constraints=ssParms)
print(ssModel$par)

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




