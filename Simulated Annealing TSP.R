# Simulated Annealing for Traveling Salesman Problem

# Function that returns the matrix of costs b/t every pair of the 17 cities
TSPCostMat <- function(){
  C <- matrix(c(0,633,257,91,412,150,80,134,259,505,353,324,70,211,268,246,121,633,0,390,661,227,
                488,572,530,555,289,282,638,567,466,420,745,518,257,390,0,228,169,112,196,154,372,
                262,110,437,191,74,53,472,142,91,661,228,0,383,120,77,105,175,476,324,240,27,182,
                239,237,84,412,227,169,383,0,267,351,309,338,196,61,421,346,243,199,528,297,150,
                488,112,120,267,0,63,34,264,360,208,329,83,105,123,364,35,80,572,196,77,351,63,0,
                29,232,444,292,297,47,150,207,332,29,134,530,154,105,309,34,29,0,249,402,250,314,
                68,108,165,349,36,259,555,372,175,338,264,232,249,0,495,352,95,189,326,383,202,236,
                505,289,262,476,196,360,444,402,495,0,154,578,439,336,240,685,390,353,282,110,324,
                61,208,292,250,352,154,0,435,287,184,140,542,238,324,638,437,240,421,329,297,314,
                95,578,435,0,254,391,448,157,301,70,567,191,27,346,83,47,68,189,439,287,254,0,145,
                202,289,55,211,466,74,182,243,105,150,108,326,336,184,391,145,0,57,426,96,268,420,
                53,239,199,123,207,165,383,240,140,448,202,57,0,483,153,246,745,472,237,528,364,
                332,349,202,685,542,157,289,426,483,0,336,121,518,142,84,297,35,29,36,236,390,238,
                301,55,96,153,336,0),nrow=17, byrow = TRUE)
  return(C)  
}

# Function that returns the total cost of a candidate solution for TSP, given the 
# proposed solution x and the cost matrix C.
totCost <- function(x, C){
  cost <- 0
  n <- length(x)
  for (i in 1:(n-1)){
    cost <- cost + C[x[i], x[i+1]]
  }
  cost <- cost + C[x[n],x[1]]
  return (cost)
}

# Function that performs Simulated Annealing optimization.
# Returns the solution vector x and the cost of the solution s
simulatedAnnealing <-function(T0, beta, iter){
  # Initialize
  T <- T0                                   # starting temperature
  C <- TSPCostMat()                         # cost matrix of distances
  n <- nrow(C)                              # number of cities
  x <- sample(1:n, size = n, replace=FALSE) # initial starting state
  Sx <- totCost(x, C)                       # cost of sequence x
  xbest <- x                                # initialize best solution vector x
  sbest <- Sx                               # initialize best solution cost Sx
  results <<- numeric(iter)                 # vector to save the cost of current solution
  
  # Perform simulated annealing procedure:
  for (i in 1:iter){
    # Choose another random permutation of the integers
    per <- sample(1:n, size = n, replace=FALSE)
    
    # Sort the first 2 in the new list
    I <- sort(per[1:2]) # vector of two indices
    
    # Create a candidate solution by reordering current x between I[1]
    # and I[2]; y is now the new proposal solution
    if (I[1]==1 & I[2]==n){
      y <- rev(x[I[1]: I[2]])
    } else if (I[1]==1 & I[2]!=n){
      y <- c(rev(x[I[1]: I[2]]), x[(I[2]+1):n])
    } else if (I[2]==n){
      y <- c(x[1:(I[1]-1)], rev(x[I[1]: I[2]]))
    } else {
      y <- c(x[1:(I[1]-1)], rev(x[I[1]: I[2]]), x[(I[2]+1):n])
    }
    
    # Compute the cost of y and determine if accepted
    Sy <- totCost(y, C)
    alpha <- ifelse(Sy < Sx, 1, exp(-(Sy-Sx)/T))
    
    # Generate ~U(0,1) and set x_t+1
    u <- runif(1)
    if (u < alpha){
      x <- y
      Sx <- Sy
    }
    
    # Select a new temperature, decrease via step-size beta
    T <- beta*T
    
    # Set new best solution vector x and cost Sx
    xbest <- x
    sbest <- Sx
    results[i] <<- sbest
  }
  return(list(solution=xbest, cost = sbest))
}


# Initialize parameters and test the code
temp0 <- 1;             
beta <- 0.9999;
numiter <- 10000;

# Get a random tour and compute its cost
set.seed(1234)
(r <- sample(1:17, size = 17, replace=FALSE))
C.tsp <- TSPCostMat()
(rCost <- totCost(r, C.tsp))

# Run simulated annealing to get the minimum solution
simulatedAnnealing(temp0, beta, numiter)

# Plot the cost of the current best solution as a function of iteration number
library(ggplot2)
ggplot(data.frame(iteration=seq(1,10000), cost=results), aes(iteration, cost)) + geom_point()+
  labs(title = "Plot of Cost of Current Best Solution vs. Iteration Number", 
       x = "Iteration Number", y = "Cost of Current Best Solution")

# Repeat the simulated annealing procedure 4 times
soln.sa <- NULL
cost.sa <- numeric(4)
instance <- c(1, 2, 3, 4)
for (i in 1:4){
  run.MCMC <- simulatedAnnealing(temp0, beta, numiter)
  cost.sa[i] <- run.MCMC$cost
  soln.sa <- rbind(soln.sa, run.MCMC$solution)
}
soln.sa <- data.frame(soln.sa)
sum.sa <- cbind(instance, soln.sa, cost.sa)
names(sum.sa) <- c("Instance", " ", " ", " ", " ", " ", " ", " ", "Tour", " ",
                   " ", " ", " ", " ", " ", " ", " ", " ", "Cost")
sum.sa
