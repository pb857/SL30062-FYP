# INITIALISE PARAMETER VARIABLES

N <- 1000 # set population size
ngen <- 8000 # set number of generations
g <- 50000 # set number of sites to index
s <- 0 # set selection coefficient for mutant
mu <- 3e-6  # set mutation rate

# Values calculated over the last X gens in a continuous range: 
record_range <- 3000:ngen
burn_range <- 1:(min(record_range) - 1)
n_intervals <- 5  # Define number of x-tick intervals (number of ticks = n_intervals + 1)

fx_tracker <- 0

# Define set of absolute fitnesses:
# abs_fit[1]: WT allele
# abs_fit[2]: mutant allele
abs_fit <- c(1, 1-s)

# make a vector with sites as rows and p(t) as columns
p_list <- rep(0, g)

# Define a FUNCTION 'mean_fitness()' to calculate mean population fitness in one generation
mean_fitness <- function(p, abs_fit) {
  q <- 1 - p  # frequency of WT allele
  w_bar <- p * abs_fit[2] + q * abs_fit[1]
  return(w_bar)
}

# Make matrices for summary parameters
px_data <- matrix(0, nrow = length(record_range), ncol = 1)  # Px
dx_data <- matrix(0, nrow = length(record_range), ncol = 1)  # Dx
h_data <- matrix(0, nrow = length(record_range), ncol = 1)  # H

# t0 GENERATION
# ====================================================================

# Allow every site to have a chance to initialise as a mutant in t0 generation

for (site in 1:g) {
  
  mutation_t0 <- rbinom(1, N, mu)
  p_list[site] <- mutation_t0 / N
  
}

# SIMULATE FROM t1 to tN
# ====================================================================

for (gen in 2:ngen){  # process for n generations
  
  fixation_counter <- 0
  
  for (site in 1:g) {   # for every site 
    
    # SELECTION IF SITE IS MUTATED
    if (p_list[site] > 0 & p_list[site] < 1){
      
      p_t <- p_list[site]
      w_bar <- mean_fitness(p_t, abs_fit)
      
      # add effect of natural selection on new frequency of p
      p_prime <- p_t * (abs_fit[2] / w_bar)
      
      # randomly sample number of mutant alleles based on p_prime
      p_next <- rbinom(1, N, p_prime) / (N)
      # set freq of mutant allele as p in current gen
      p_list[site] <- p_next
      
    }
    
    # RETURN FIXED SITES TO MUTABLE POOL 
    else if (p_list[site] == 1){
      fixation_counter <- fixation_counter + 1
      p_list[site] <- 0
      
    }
    
    # MUTATION IF SITE IS UNMUTATED
    else if (p_list[site] == 0){
      mutation_event <- rbinom(1, N, mu) / (N)
      p_list[site] <- mutation_event
      
    }
    
  } # end of inner loop
  
  # Passively record fixations during burn-in period
  if (gen %in% burn_range){
    fx_tracker <- fx_tracker + fixation_counter
  }
  
  # Active recording of representative samples of summary parameters in stationary distribution
  if (gen %in% record_range){
    record_index <- gen - min(record_range) + 1
    
    px_data[record_index, 1] <- sum(p_list > 0 & p_list < 1) / (g)
    
    h_temp <- p_list * (1 - p_list)
    h_data[record_index, 1] <- mean(h_temp)
    
    if (gen == min(record_range)){
      dx_data[record_index, 1] <- (fx_tracker + fixation_counter) / g
    }
    else{
      dx_data[record_index, 1] <- dx_data[record_index-1, 1] + (fixation_counter / g)
    }
    
  }
  
  # Print progress every 500 generations
  if (gen %% 500 == 0) {
    cat("Completed generation:", gen, "\n")
  }
  
}

# ==================================================================

# Expected number of segregating sites
theta <- 2*N*mu
sum <- sum(1/(1:(N-1))) 
E_S <- theta * sum

# Correction for E[S]
gE <- g - (E_S * g)   # effective number of mutable sites 
muE <- mu * (gE / g)  # effective mutation rate 
thetaE <- 2*N*muE     # effective population-scaled mutation rate
SE <- thetaE * sum    # E[S] corrected

# Expected corrected level of heterozygosity
E_H <- thetaE / (1 + thetaE)

# PLOTS 
# =====

# Dynamically calculate x-ticks based on record_range
x_ticks <- seq(from = min(record_range), 
               to = max(record_range), 
               length.out = n_intervals + 1)  # +1 to include both ends

# PX:

matplot(record_range, px_data, type = "l", col = "red",
        xaxt = "n",   # Suppress x-axis
        xlab = "Generations", 
        ylab = "Px",
        main = "Px over time")

axis(1, at = x_ticks, labels = round(x_ticks))  

abline(h = E_S , col = "black", lwd = 1.5, lty = 2) # draw raw E[S] as a reference line 
abline(h = SE , col = "green", lwd = 1.5, lty = 2) # draw corrected E[S] as a reference line 
abline(h = mean(px_data), col = "blue", lwd = 1.5, lty = 2) # empirical mean Px

legend("bottomleft", 
       legend = c(paste("N:", N),
                  paste("ngen", ngen),
                  paste("g:", g),
                  paste("mu:", mu)
       ),
       bty = "n"
)

legend("bottomright", 
       legend = c(paste("Raw E[S]:", E_S),
                  paste("Corrected E[S]:", SE),
                  paste("Mean Px:", mean(px_data))
       ),
       bty = "n"
)

# DX: 

matplot(record_range, dx_data, type = "l", col = "blue",
        xaxt = "n",   # Suppress x-axis
        xlab = "Generations", 
        ylab = "Dx",
        main = "Dx over time")

axis(1, at = x_ticks, labels = round(x_ticks))  

legend("bottomright", 
       legend = c(paste("N:", N),
                  paste("ngen", ngen),
                  paste("g:", g),
                  paste("mu:", mu),
                  paste("Raw E[S]:", E_S),
                  paste("Corrected E[S]:", SE)
       ),
       bty = "n"
)

# H:

matplot(record_range, h_data, type = "l", col = "green",
        xaxt = "n",
        xlab = "Generations", 
        ylab = "H",
        main = "H over time ")

axis(1, at = x_ticks, labels = round(x_ticks)) 

legend("bottomright", 
       legend = c(paste("N:", N),
                  paste("ngen", ngen),
                  paste("g:", g),
                  paste("mu:", mu),
                  paste("Corrected E[H]:", E_H)
       ),
       bty = "n"
)

# =================================================================