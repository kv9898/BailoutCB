######## INFO ########

# PROJECT
# Replication files for: Bail out the money printer? The impact of fiscal indemnity on central bank profitability

# Authors: Dianyi Yang
# R Script
# Purpose: This script generates the numerical evaluation results for the formal model.
# Created: 24 JUL 2024
# Inputs: None
# Outputs: stored_output/formal_results.rds

######## SETUP ########
#setwd("~/BailoutCB") # change the working directory to your folder

need <- c('tidyverse') # list packages needed
have <- need %in% rownames(installed.packages()) # checks packages you have
if(any(!have)) install.packages(need[!have]) # install missing packages
invisible(lapply(need, library, character.only=T))

# initialise parameters
rho <- seq(1.01,1.99, length.out=396) #1.376
alpha <- seq(0.01,4*sqrt(3), length.out=400) #4*sqrt(3)
psi=1
results <- data.frame(rho=rep(rho, each=length(alpha)), alpha=rep(alpha, length(rho)))
results <- results %>% mutate(
  r1 = (-2 * psi - 1 + sqrt((2 * psi + 1)^2 - 4 * psi * ( 1 - rho ))) / (2 * psi)
)
ggplot(results, aes(x=rho, y=r1)) + geom_line() # inspect the relationship between r(x=1) and rho

# set bound
results <- results |>
  mutate(
    within = ifelse(alpha <= (4 - 2 * rho) / (psi * r1), T, F)
  )

results <- results |>
  mutate(
    r1 = ifelse(within, r1, NA)
  )

ggplot(results, aes(x=rho, y=alpha, fill=within)) + geom_tile() # check the shape of bound

# optimise over r0
get_r0 <- function(row) {
 if (row["within"] == F) {
   return(NA)
 }
 rho <- row["rho"] |> as.numeric()
 alpha <- row["alpha"] |> as.numeric()
 r1 <- row["r1"] |> as.numeric()
 gradient_0 = ( - (2 - rho)^2 / 8 ) * alpha * psi - 2 * (1- rho)
 constraint = (4 - 2 * rho) / (alpha * psi)
 if (gradient_0<0){
   r0 = ifelse(r1<=constraint, 0, 99999)
 } else{
   a= (alpha * psi )^3 / 16 + 4 * psi
   b= (2 - rho) / 8 * (alpha * psi)^2 - 4 * psi - 2
   max = (b + sqrt( b^2 - 2 * a * ( (2-rho)^2 / 8 * alpha * psi + 2 - 2 * rho ) )) / a
   if ( r1 > constraint) { #max > constraint) {
     r0 = r1
   } else{
     r0 = min(r1, max)
   }
 }
 return(r0)
}
results$r0 <- apply(results, 1, get_r0)
ggplot(results, aes(x=rho, y=alpha, fill=(r0<r1))) + geom_tile() # r(x=0)<r(x=1) for our domain


#calculate expected utility of government with indemnity
  #recall that inflation is always 2 and expected growth is also 2
  #recall that the government's utility is only determined by growth in equilibrium
#EUG1=4
#ggplot(results, aes(x=rho, y=alpha, fill=EUG1)) + geom_tile() + scale_fill_viridis_c()

#calculate expected utility of government without indemnity
get_EUG0 <- function(row) {
  if (row["within"] == F) {
    return(NA)
  }
  rho <- row["rho"] |> as.numeric()
  alpha <- row["alpha"] |> as.numeric()
  r0 <- row["r0"] |> as.numeric()
  UG1 = 3-rho+r0+psi * (r0 + 2) * r0
  if ((alpha * psi * r0)/2 > 2-rho ) {
    EUG2 = 1 / 8 * rho^2 - 1 / 2 * rho + 5 / 2
  } else {
    EUG2 = - (alpha * psi * r0)^2 / 32 + (alpha * psi * r0) / 4 - (alpha * psi * r0) * rho / 8 + 2
  }
  return (UG1+EUG2)
}

results$EUG0 <- apply(results, 1, get_EUG0) # calculate EU_G(x=0)

#####################Period 2 evalulation###############
results <- results |> mutate(
  i1 = (4 - rho + 2 * psi * r1 ) / (psi * r1 + 1),
  i0 = case_when(
    4 - rho <= 2 + alpha*psi*r0/2 ~ 2,
    4 - rho > 2 + alpha*psi*r0/2 ~ 4 - rho - (alpha * psi * r0) /2
  ),
  p1 = -psi * (i1-2) * r1,
  p0 = -psi * (i0-2) * r0,
  p1_pre = -psi * (-r1-2) * r1,,
  p0_pre = -psi * (-r0-2) * r0,
)

results <- results |> mutate(
  pi0=6-rho-i0
)
# save results
saveRDS(results, "stored_output/formal_results.rds")