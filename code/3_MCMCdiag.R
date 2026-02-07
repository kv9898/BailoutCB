######## INFO ########

# PROJECT
# Replication files for: Bail out the money printer? The impact of fiscal indemnity on central bank profitability

# Authors: Dianyi Yang
# R Script
# Purpose: This script is responsible for the MCMC diagnosis of Main Analysis.
# Created: 13 Aug 2024
# Inputs: stored_output/result.profit.full.RData
# Outputs: stored_output/MCMC.RData

######## SETUP ########
#setwd("~/BailoutCB") # change the working directory to your folder

library(coda)
library(kableExtra)

load("stored_output/result.profit.full.RData")

mc_pre <- sapply(1:10000, function(i) {mean(out.full.cov$yo_t[11:23] - out.full.cov$yct[11:23, i])}) |> mcmc()
mc_post <- sapply(1:10000, function(i) {mean(out.full.cov$yo_t[24:25] - out.full.cov$yct[24:25, i])})|> mcmc()

#summary statistics
sum <- rbind(t(summary(mc_pre)$statistics), t(summary(mc_post)$statistics))
sum <- cbind('Period'=c('Deflationary','Inflationary'),
              "Start year" = c("2009","2022"),
              "End year" = c("2021","2023"),
              as.data.frame(sum))
sum_mcmc <- kbl(sum, centering=T, digits=6, format="latex") |> kable_styling(latex_options="scale_down")


#quantiles
quant <- rbind(t(summary(mc_pre)$quantiles), t(summary(mc_post)$quantiles))
quant <- cbind('Period'=c('Deflationary','Inflationary'),
              "Start year" = c("2009","2022"),
              "End year" = c("2021","2023"),
              as.data.frame(quant))

sum_quant <- kbl(quant, centering=T, digits=6, format="latex") |> kable_styling(latex_options="scale_down")

save(mc_pre, mc_post, sum_mcmc, sum_quant, file='stored_output/MCMC.RData')
