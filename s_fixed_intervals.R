# INITIALISE PARAMETER VARIABLES

# install.packages("writexl")
library(writexl)

N <- 1000  # set population size 
ngen <- 28000  # set number of generations
g <- 50000   # set number of sites to index
s_range <- seq(0, 0.005, by = 0.0001)  # set range and interval of selection coefficients for mutant
mu <- 3e-6  # set mutation rate

# Values calculated over the last X gens in a continuous range: 
record_range <- 3000:ngen # manually change starting index to change recording range
burn_range <- 1:(min(record_range) - 1) # burn-in period for which number of fixations is passively counted up

# ===================================================================================================================

# Define a FUNCTION 'mean_fitness()' to calculate mean population fitness in one generation
mean_fitness <- function(p, abs_fit) {
  q <- 1 - p  # frequency of WT allele
  w_bar <- p * abs_fit[2] + q * abs_fit[1]
  return(w_bar)
}

# create a matrix to store average Px, Dx and H for each selection coefficient in s_range

s_data <- matrix(0, nrow = length(s_range), ncol = 3)
rownames(s_data) <- as.character(s_range)
colnames(s_data) <- c("Px", "Dx", "H")

# MAIN LOOP
# ===================================================================================================================

for (s in s_range) {
  
  # make a vector with sites as rows and p(t) as columns
  p_list <- rep(0, g)
  
  abs_fit <- c(1, 1-s)  # initialise set of absolute fitness for non-mutant and mutant 
  fx_tracker <- 0 # initialise number of fixations as 0 in every run of s (passive burn-in period tracking)
  
  # Initialise summary parameter lists for every iteration of s
  # populate column matrices with ngen rows of zeroes
  px_data <- matrix(0, nrow = length(record_range), ncol = 1)  # Px
  dx_data <- matrix(0, nrow = length(record_range), ncol = 1)  # Dx
  h_data <- matrix(0, nrow = length(record_range), ncol = 1)  # H
  
  # Allow every site to have a chance to initialize as a mutant in t0 generation
  for (site in 1:g) {
    mutation_t0 <- rbinom(1, N, mu)
    p_list[site] <- mutation_t0 / N
  }
  
  # Simulate generations t1 -> tN
  for (gen in 2:ngen) {
    
    fixation_counter <- 0 # tracks new number of fixations per generation
    
    for (site in 1:g) {
      
      # Have site undergo SELECTION if site is mutated (0 < site < 1)
      if (p_list[site] > 0 & p_list[site] < 1) {
        p_t <- p_list[site]
        w_bar <- mean_fitness(p_t, abs_fit)
        
        p_prime <- p_t * (abs_fit[2] / w_bar) # add effect of selection to new p frequency
        p_next <- rbinom(1, N, p_prime) / N   # randomly sample number of alleles based on p_prime
        p_list[site] <- p_next                # set frequency of mutant allele as p in current gen
      } 
      
      # Return any fixed sites (site == 1) to the mutable pool by resetting their frequency to zero
      else if (p_list[site] == 1) {
        fixation_counter <- fixation_counter + 1  # count up number of fixations by 1
        p_list[site] <- 0 # reset site to zero to make it mutable
      } 
      
      # Chance for un-mutated site (site == 0) to undergo MUTATION
      else if (p_list[site] == 0) {
        mutation_event <- rbinom(1, N, mu) / N
        p_list[site] <- mutation_event
      }
      
    } # End of site loop 
    
    # Passively record fixations during burn-in period
    if (gen %in% burn_range){
      fx_tracker <- fx_tracker + fixation_counter
    }
    
    # Active recording of representative samples of summary parameters in stationary distribution
    if (gen %in% record_range){
      record_index <- gen - min(record_range) + 1
      
      px_data[record_index, 1] <- sum(p_list > 0) / (g)
      
      h_temp <- p_list * (1 - p_list)
      h_data[record_index, 1] <- mean(h_temp)
      
      if (gen == min(record_range)){
        dx_data[record_index, 1] <- (fx_tracker + fixation_counter) / g
      }
      else{
        dx_data[record_index, 1] <- dx_data[record_index-1, 1] + (fixation_counter / g)
      }
      
    }
    
  } # End of generation loop
  
  # storage of representative samples of each summary parameter Px, Dx and H
  s_data[as.character(s), 1] <- mean(px_data)
  s_data[as.character(s), 2] <- mean(dx_data)
  s_data[as.character(s), 3] <- mean(h_data)
  
  cat("Completed s =", s, "at index =", which(s_range == s), "\n\n")
  
}

# ==================================================================

# EXPECTED PARAMETERS AT NEUTRALITY
# =================================

# Expected number of segregating sites
theta <- 2*N*mu
sum <- sum(1/(1:(N-1))) 
E_S <- theta * sum

# Correction for E[S]
gE <- g - (E_S*g) # number of mutable sites 
muE <- mu * (gE / g)  # effective mutation rate 

thetaE <- 2*N*muE
SE <- thetaE * sum

# Equilibrium level of neutral heterozygosity
E_H <- thetaE / (1 + thetaE)

# PLOTS 
# ==================================================================

#par(mfrow = c(1,2)) # panelise the plots onto one window

# Px: 
matplot(1:nrow(s_data), s_data[, 1], type = "l", col = "red", xaxt = "n",
        xlab = "s", 
        ylab = "Px",
        main = "Px-s ")
axis(1, at = 1:nrow(s_data), labels = rownames(s_data))
grid()
legend("topright",
       legend = c(paste("N:", N),
                  paste("g:", g),
                  paste("gens:", ngen),
                  paste("mu:", mu),
                  paste("E[S]:", SE),
                  paste("Recorded over last", ngen - min(record_range), "gens")),
       bty = "n"
)

# Dx: 
matplot(1:nrow(s_data), s_data[, 2], type = "l", col = "blue", xaxt = "n",
        xlab = "s", 
        ylab = "Dx",
        main = "Dx-s ")
axis(1, at = 1:nrow(s_data), labels = rownames(s_data))
grid()
legend("topright",
       legend = c(paste("N:", N),
                  paste("g:", g),
                  paste("gens:", ngen),
                  paste("mu:", mu),
                  paste("E[S]:", SE),
                  paste("Recorded over last", ngen - min(record_range), "gens")),
       bty = "n"
)

# H:
matplot(1:nrow(s_data), s_data[, 3], type = "l", col = "green", xaxt = "n",
        xlab = "s", 
        ylab = "H",
        main = "H-s ")
axis(1, at = 1:nrow(s_data), labels = rownames(s_data))
grid()
legend("topright",
       legend = c(paste("N:", N),
                  paste("g:", g),
                  paste("gens:", ngen),
                  paste("mu:", mu),
                  paste("E[H]:", E_H),
                  paste("Recorded over last", ngen - min(record_range), "gens")),
       bty = "n"
)

# =====================================================================================
# save data into an Excel file

dos_data <- as.data.frame(s_data)
write_xlsx(dos_data,"fixed_interval_run.xlsx")

# =====================================================================================
