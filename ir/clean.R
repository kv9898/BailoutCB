######## INFO ########

# PROJECT
# Replication files for: Bail out the money printer? The impact of fiscal indemnity on central bank profitability

# Authors: Dianyi Yang
# R Script
# Purpose: This script cleans the raw interest rate (ir) data, and stores the 
#          cleaned ir data in the .csv format. 
#          The cleaned data are then manually added to "cleaned_data.csv"
# Created: 04 Feb 2024
# Updated: 02 Sep 2024
# Packages used: tidyverse (2.0.0)
# Inputs: ir/raw/BoE.csv
#         ir/raw/ECB.csv
#         ir/raw/BoJ.csv
#         ir/raw/Fed.csv
#         ir/raw/SNB.csv
#         ir/raw/SWE.csv
#         ir/raw/ISR.csv

# Outputs: ir/cleaned/BoE_cleaned.csv
#          ir/cleaned/ECB_cleaned.csv
#          ir/cleaned/BoJ_cleaned.csv
#          ir/cleaned/fed_cleaned.csv
#          ir/cleaned/SNB_cleaned.csv
#          ir/cleaned/SWE_cleaned.csv
#          ir/cleaned/ISR_cleaned.csv

######## SETUP ########
#setwd("~/BailoutCB") # uncomment and change the working directory to your folder

need <- c('tidyverse') # list packages needed
have <- need %in% rownames(installed.packages()) # checks packages you have
if(any(!have)) install.packages(need[!have]) # install missing packages
invisible(lapply(need, library, character.only=T))

Sys.setlocale('LC_TIME', 'C') #set time format

######## CLEAN ########
#boe
boe_raw <- read_csv('ir/raw/BoE.csv')
boe_raw$policy_change <- as.Date(boe_raw$policy_change, format='%d %b %y')

# Assuming the rates change at the end of the policy_change date, we calculate the days in effect.
boe_raw$days_in_effect <- c(diff(-as.numeric(boe_raw$policy_change))-1, NA)

# function for calculating the weighted average interest rate for given year
weighted_average <- function(year){
  rates <- subset(boe_raw[1:which(boe_raw$policy_change<as.Date(paste0(year,'-01-01')))[1],], policy_change < as.Date(paste0(year+1,'-01-01')))
  if(nrow(rates) == 1) {
    return(rates$Rate)
  } else {
  rates$days_in_effect[1] <- as.numeric(as.Date(paste0(year+1,'-01-01')) - rates$policy_change[1])
  rates$days_in_effect[nrow(rates)] <- -as.numeric( as.Date(paste0(year,'-01-01')) - rates$policy_change[nrow(rates)-1] )
  
  weighted_sum <- sum(rates$Rate * rates$days_in_effect)
  total_days <- sum(rates$days_in_effect)
  average_rate <- weighted_sum / total_days
  return(average_rate)
  }
}

years_of_interest <- 1996:2023
boe_cleaned <- data.frame(year = years_of_interest, avg_rate = map_dbl(years_of_interest, weighted_average))
save(boe_cleaned, file='ir/cleaned/boe_cleaned.RData') # save cleaned data
write_csv(boe_cleaned, 'ir/cleaned/boe_cleaned.csv')

# ECB
ecb_raw <- read_csv('ir/raw/ECB.csv')
years_of_interest <- 1999:2023
ecb_avg <- function(year){ # only need to calculate the averages as the data is daily
  subset <- subset(ecb_raw, year(DATE) == year)
  return(mean(subset$ECBDFR))
}
ecb_cleaned <- data.frame(year = years_of_interest, avg_rate = map_dbl(years_of_interest, ecb_avg))
save(ecb_cleaned, file='ir/cleaned/ecb_cleaned.RData') # save cleaned data
write_csv(ecb_cleaned, 'ir/cleaned/ecb_cleaned.csv')

# BoJ
boj_raw <- read_csv('ir/raw/BoJ.csv')
years_of_interest <- 1996:2023
boj_avg <- function(year){ # only need to calculate the averages as the data is monthly
  subset <- subset(boj_raw, year(DATE) == year)
  return(mean(subset$RATE))
}
boj_cleaned <- data.frame(year = years_of_interest, avg_rate = map_dbl(years_of_interest, boj_avg))
write_csv(boj_cleaned, 'ir/cleaned/boj_cleaned.csv') # save cleaned data

# Fed
fed_raw <- read_csv('ir/raw/Fed.csv')
years_of_interest <- 1996:2023
fed_avg <- function(year){ # only need to calculate the averages as the data is monthly
  subset <- subset(fed_raw, year(DATE) == year)
  return(mean(subset$FEDFUNDS))
}
fed_cleaned <- data.frame(year = years_of_interest, avg_rate = map_dbl(years_of_interest, fed_avg))
write_csv(fed_cleaned, 'ir/cleaned/fed_cleaned.csv') # save cleaned data

#SNB
SNB_raw <- read_csv('ir/raw/SNB.csv')
years_of_interest <- 1998:2023
SNB_avg <- function(year){ # only need to calculate the averages as the data is monthly
  subset <- subset(SNB_raw, year(DATE) == year)
  return(mean(subset$RATE))
}
SNB_cleaned <- data.frame(year = years_of_interest, avg_rate = map_dbl(years_of_interest, SNB_avg))
write_csv(SNB_cleaned, 'ir/cleaned/SNB_cleaned.csv') # save cleaned data

#SWE
SWE_raw <- read_csv('ir/raw/SWE.csv')
SWE_raw$date <- c(as.Date(SWE_raw$DATE[1:780], format = "%Y/%m/%d"), 
      as.Date(SWE_raw$DATE[781:nrow(SWE_raw)], format = "%d/%m/%Y") ) 
      # data before and after 2020 are in different formats
      # reformat as Date object
years_of_interest <- 1999:2023
SWE_avg <- function(year){ # only need to calculate the averages as the data is monthly
  subset <- subset(SWE_raw, year(date) == year)
  return(mean(subset$RATE))
}
SWE_cleaned <- data.frame(year = years_of_interest, avg_rate = map_dbl(years_of_interest, SWE_avg))
write_csv(SWE_cleaned, 'ir/cleaned/SWE_cleaned.csv') # save cleaned data

#ISR
ISR_raw <- read_csv('ir/raw/ISR.csv')
years_of_interest <- 1999:2023
ISR_avg <- function(year){ # only need to calculate the averages as the data is tri-monthly
  subset <- subset(ISR_raw, year(DATE) == year)
  return(mean(subset$RATE))
}
ISR_cleaned <- data.frame(year = years_of_interest, avg_rate = map_dbl(years_of_interest, ISR_avg))
write_csv(ISR_cleaned, 'ir/cleaned/ISR_cleaned.csv') # save cleaned data
