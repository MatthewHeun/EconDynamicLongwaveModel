# Load raw data
rawHistoricalData <- read.table(file="data/USData.txt", header=TRUE, row.names=1, sep="\t")
# Duplicate the raw data into a new data.frame that we'll extend with calculated historical data.
calculatedHistoricalData <- rawHistoricalData
# Calculate other raw data items
L_Y <- calculatedHistoricalData["L"] - calculatedHistoricalData["L_A"]; names(L_Y)  <- "L_Y"
calculatedHistoricalData <- cbind(calculatedHistoricalData, L_Y)
# Calculate indexed data
firstYear <- "1980"
w_worker_0 <- calculatedHistoricalData[firstYear, "w_worker"] # Initial wage
Y_0  <- calculatedHistoricalData[firstYear, "Y"]  # Initial GDP
N_0 <- calculatedHistoricalData[firstYear, "N"]   # Initial population
L_Y_0 <- L_Y[firstYear,"L_Y"] # Initial L_Y
y <- calculatedHistoricalData["Y"] / Y_0; names(y) <- "y" # Indexed GDP
n <- calculatedHistoricalData["N"] / N_0; names(n) <- "n" # Indexed population
l_Y <- (calculatedHistoricalData["L_Y"]) / L_Y_0; names(l_Y) <- "l_Y" # Indexed workers in the manufacturing sector.
calculatedHistoricalData <- cbind(calculatedHistoricalData, y, n, l_Y)

# Calculate per-capita information
# Births/year-capita
b <- calculatedHistoricalData["B"] / calculatedHistoricalData["N"]
names(b) <- "b" 
b_0 <- b[firstYear, "b"]
# Deaths/year-capita
d <- calculatedHistoricalData["D"] / calculatedHistoricalData["N"]
names(d) <- "d" 
d_0 <- d[firstYear, "d"]
calculatedHistoricalData <- cbind(calculatedHistoricalData, b, d)

# Fraction of workers in each year
tau <- calculatedHistoricalData["L"] / calculatedHistoricalData["N"]; 
names(tau) <- "tau" 
# consumption in 2005$/year-capita. Same as wages in 2005$/year-capita.
c <- calculatedHistoricalData["w_worker"] * tau; 
names(c) <- "c" 
calculatedHistoricalData <- cbind(calculatedHistoricalData, tau, c)

# Calculate nu, pi, and zeta
# **** SEE JonesModelGroundTruthParameters.EES IN THE REPOSITORY ****
nu <- calculatedHistoricalData["b"] / (1.0 - calculatedHistoricalData["tau"]); 
names(nu) <- "nu" # Historical nu
nu_0 <- nu[firstYear, "nu"]
pi <- calculatedHistoricalData["L_A"] / calculatedHistoricalData["L"]
names(pi) <- "pi" # Historical pi
pi_0 <- pi[firstYear, "pi"]
wALA <- calculatedHistoricalData["w_worker"] * calculatedHistoricalData["L_A"]; 
names(wALA) <- "wALA"
piY <- pi * calculatedHistoricalData["Y"]; names(piY) <- "piY"
zeta <- wALA / piY; names(zeta) <- "zeta"
zeta_0 <- zeta[firstYear, "zeta"]
calculatedHistoricalData <- cbind(calculatedHistoricalData, nu, pi, zeta)

# calculate b_bar, b_tilde, c_bar, c_tilde, d_bar, z
# Asymptotic birth rate as consumption goes to infinity
b_bar = 0.0
# Birth rate above asymptote, b_tilde
b_tilde <- calculatedHistoricalData["b"] - b_bar
names(b_tilde) <- "b_tilde"
b_tilde_0 <- b_tilde[firstYear, "b_tilde"]
# Consumption for subsistence level
############## Perform this calculation here #################
c_bar = 5655 # 2005$/year
# Consumption above subsistence level, c_tilde
c_tilde <- calculatedHistoricalData["c"] - c_bar
names(c_tilde) <- "c_tilde"
c_tilde_0 <- c_tilde[firstYear, "c_tilde"]
# Asymptotic mortality rate
d_bar <- 0
# calculate z
z <- calculatedHistoricalData["c"] / c_bar - 1
names(z) <- "z" 
z_0 <- z[firstYear, "z"]
calculatedHistoricalData <- cbind(calculatedHistoricalData, b_tilde, c_tilde, z)
