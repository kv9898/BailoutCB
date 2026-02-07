######## INFO ########

# PROJECT
# Replication files for: Bail out the money printer? The impact of fiscal indemnity on central bank profitability

# Authors: Dianyi Yang
# R Script
# Purpose: This script is responsible for the empirical analysis of the paper,
#          including, that in the appendices, except for the MCMC diagnosis.
# Created: 09 Aug 2024
# Inputs: cleaned_data.rds
#         custom_functions/DM-LFM.r
#         custom_functions/solver.R
# Outputs: stored_output/profit.placebo.RData
#          stored_output/result.profit.full.RData
#          stored_output/profit.prepost.RData
#          stored_output/profit.early.RData
#          stored_output/result.profit.early.RData
#          stored_output/supp-ir.RData
#          stored_output/supp-bs.RData
#          stored_output/suppresults.RData
#          stored_output/alt.RData

# SETUP --------------------------------------------------------------------
# setwd("~/BailoutCB") # change the working directory to your folder

library(bpCausal) # https://github.com/liulch/bpCausal
library(tidyverse)
library(modelsummary)
library(tinytable)

# READ AND MANIPULATE DATA -------------------------------------------------
empdata <- read_rds("cleaned_data.rds")
source("custom_functions/DM-LFM.r")

# TEST OF PRETREND ------------------------------------------------------------
data_pre <- empdata |> filter(year <= 2008)
data_pre <- data_pre |>
  mutate(
    treat = if_else(Bank == "BoE" & year >= 2004, 1, 0),
  )

set.seed(666)

data_pre <- data_pre |>
  group_by(Bank) |>
  arrange(year, .by_group = TRUE) |>
  ungroup()

class(data_pre) <- "data.frame"

out.pre <- bpCausal(
  data = data_pre,
  Yname = "profit_ppgdp",
  Dname = "treat",
  Xname = c(),
  # Zname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Zname = c("traded", "Eurozone", "reappointable"),
  # Aname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Aname = c("traded", "Eurozone", "reappointable"),
  index = c("Bank", "year"),
  re = "both",
  ar1 = TRUE,
  r = 10,
  niter = 15000,
  burn = 5000,
  xlasso = 1,
  zlasso = 1,
  alasso = 1,
  flasso = 1,
  a1 = 0.001,
  a2 = 0.001,
  b1 = 0.001,
  b2 = 0.001,
  c1 = 0.001,
  c2 = 0.001,
  p1 = 0.001,
  p2 = 0.001
) # DM-LFM Warning: a bit time consuming
out.pre.nocov <- bpCausal( # without covariates
  data = data_pre,
  Yname = "profit_ppgdp",
  Dname = "treat",
  Xname = c(),
  # Zname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Zname = c(),
  # Aname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Aname = c(),
  index = c("Bank", "year"),
  re = "both",
  ar1 = TRUE,
  r = 10,
  niter = 15000,
  burn = 5000,
  xlasso = 1,
  zlasso = 1,
  alasso = 1,
  flasso = 1,
  a1 = 0.001,
  a2 = 0.001,
  b1 = 0.001,
  b2 = 0.001,
  c1 = 0.001,
  c2 = 0.001,
  p1 = 0.001,
  p2 = 0.001
) # DM-LFM Warning: a bit time consuming

result.use.pre <- effSummary(out.pre) # extract effect information
placebo.plot <- ggplot(result.use.pre$est.eff, aes(x = time + 2003)) +
  geom_ribbon(aes(ymin = counterfactual_ci_l, ymax = counterfactual_ci_u, fill = "95% CI"), alpha = 0.5) +
  geom_line(aes(y = observed, color = "Observed"), linewidth = 0.75, alpha = 0.5) +
  geom_line(aes(y = estimated_counterfactual, color = "Counterfactual"), linewidth = 0.75, alpha = 0.5) +
  geom_vline(xintercept = 2003, linetype = "dashed") +
  geom_vline(xintercept = 2008) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_color_manual(
    name = "",
    breaks = c("Observed", "Counterfactual"),
    values = c("Observed" = "#00BFC4", "Counterfactual" = "#F8766D")
  ) +
  labs(x = "Year", y = "BoE profit(% GDP)") +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(1999, 2008, 1), minor_breaks = NULL) +
  ylim(-0.5, 0.5) +
  scale_fill_manual(
    name = c("", ""),
    breaks = c("95% CI"),
    values = c("lightgrey")
  )
profit.placebo <- coefCustom(out.pre, custom.period = 1:4)
profit.placebo[1, ] <- effCustom(out.pre, custom.period = 1:4)[1, ]
rownames(profit.placebo) <- c("ATT", "traded", "Eurozone", "reappointable")

profit.placebo <- list(results = profit.placebo, obs = nrow(data_pre), tr = 1, ct = length(unique(data_pre$Country)) - 1)
class(profit.placebo) <- "DMLFM"

profit.placebo.nocov <- effCustom(out.pre.nocov, custom.period = 1:4)
rownames(profit.placebo.nocov) <- c("ATT")
profit.placebo.nocov <- list(results = profit.placebo.nocov, obs = nrow(data_pre), tr = 1, ct = length(unique(data_pre$Country)) - 1)
class(profit.placebo.nocov) <- "DMLFM"


tidy.DMLFM <- function(x, ...) {
  ret <- data.frame(
    term      = x$results |> rownames(),
    estimate  = x$results[, "mean"],
    conf.low  = x$results[, "ci_l"],
    conf.high = x$results[, "ci_u"],
    p.value   = x$results[, "p"]
  )
  ret
}
glance.DMLFM <- function(x, ...) {
  ret <- data.frame(
    nobs = x$obs,
    tr = x$tr,
    ct = x$ct
  )
  ret
}

gof_map <- list(
  list("raw" = "nobs", "clean" = "Observations", "fmt" = 0),
  list("raw" = "tr", "clean" = "Treated Units", "fmt" = 0),
  list("raw" = "ct", "clean" = "Control Units", "fmt" = 0)
)
cm <- c(
  "ATT" = "ATT",
  "traded" = "Publicly traded",
  "Eurozone" = "Euro Area",
  "reappointable" = "Reappointability"
)

placebo.table <- modelsummary(list(profit.placebo.nocov, profit.placebo),
  output = "latex",
  estimate = "{estimate}{stars}",
  statistic = "[{conf.low}, {conf.high}]",
  gof_map = gof_map, coef_map = cm, escape = F,
  notes = c(
    "+ p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001",
    "95\\% equal-tailed Credible Intervals in square brackets."
  )
)
save(profit.placebo, placebo.plot, placebo.table, file = "stored_output/profit.placebo.RData")

# full data -----------------------------------------------------------------
## with covariates --------------------------------------------------------
set.seed(123)
empdata <- empdata |>
  mutate(
    treat = if_else(Bank == "BoE" & year >= 2009, 1, 0),
  ) |>
  filter(!Country %in% c("Switzerland", "Eurozone"))
class(empdata) <- "data.frame"
out.full.cov <- bpCausal(
  data = empdata,
  Yname = "profit_ppgdp",
  Dname = "treat",
  Xname = c(),
  # Zname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Zname = c("traded", "Eurozone", "reappointable"),
  # Aname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Aname = c("traded", "Eurozone", "reappointable"),
  index = c("Bank", "year"),
  re = "both",
  ar1 = TRUE,
  r = 10,
  niter = 15000,
  burn = 5000,
  xlasso = 1,
  zlasso = 1,
  alasso = 1,
  flasso = 1,
  a1 = 0.001,
  a2 = 0.001,
  b1 = 0.001,
  b2 = 0.001,
  c1 = 0.001,
  c2 = 0.001,
  p1 = 0.001,
  p2 = 0.001
) # DM-LFM Warning: a bit time consuming

result.use.full.cov <- effSummary(out.full.cov) # extract effect information
save(out.full.cov, result.use.full.cov, file = "stored_output/result.profit.full.RData")

profit.pre.cov <- coefCustom(out.full.cov, custom.period = 1:13)
profit.pre.cov[1, ] <- effCustom(out.full.cov, custom.period = 1:13)[1, ]
rownames(profit.pre.cov) <- c("ATT", "traded", "Eurozone", "reappointable")

profit.post.cov <- coefCustom(out.full.cov, custom.period = 14:15)
profit.post.cov[1, ] <- effCustom(out.full.cov, custom.period = 14:15)[1, ]
rownames(profit.post.cov) <- c("ATT", "traded", "Eurozone", "reappointable")

profit.pre.cov <- list(results = profit.pre.cov, obs = nrow(empdata), tr = 1, ct = length(unique(empdata$Country)) - 1)
profit.post.cov <- list(results = profit.post.cov, obs = nrow(empdata), tr = 1, ct = length(unique(empdata$Country)) - 1)
class(profit.pre.cov) <- "DMLFM"
class(profit.post.cov) <- "DMLFM"

profit.2023 <- effCustom(out.full.cov, custom.period = 15)[1, ] # 2023 results are mentioned in the manuscript

## without covariates -----------------------------------------------------
set.seed(123)
class(empdata) <- "data.frame"
out.full <- bpCausal(
  data = empdata,
  Yname = "profit_ppgdp",
  Dname = "treat",
  Xname = c(),
  # Zname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Zname = c(),
  # Aname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Aname = c(),
  index = c("Bank", "year"),
  re = "both",
  ar1 = TRUE,
  r = 10,
  niter = 15000,
  burn = 5000,
  xlasso = 1,
  zlasso = 1,
  alasso = 1,
  flasso = 1,
  a1 = 0.001,
  a2 = 0.001,
  b1 = 0.001,
  b2 = 0.001,
  c1 = 0.001,
  c2 = 0.001,
  p1 = 0.001,
  p2 = 0.001
) # DM-LFM Warning: a bit time consuming
profit.pre <- effCustom(out.full, custom.period = 1:13)
rownames(profit.pre) <- c("ATT")
profit.post <- effCustom(out.full, custom.period = 14:15)
rownames(profit.post) <- c("ATT")

profit.pre <- list(results = profit.pre, obs = nrow(empdata), tr = 1, ct = length(unique(empdata$Country)) - 1)
profit.post <- list(results = profit.post, obs = nrow(empdata), tr = 1, ct = length(unique(empdata$Country)) - 1)
class(profit.pre) <- "DMLFM"
class(profit.post) <- "DMLFM"


## main table ------------------------------------------------------------
profit.table <- modelsummary(list(profit.pre, profit.pre.cov, profit.post, profit.post.cov),
  estimate = "{estimate}{stars}", output = "latex",
  statistic = "[{conf.low}, {conf.high}]",
  gof_map = gof_map, coef_map = cm, escape = F,
  notes = c(
    "+ p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001",
    "95\\% equal-tailed Credible Intervals in square brackets."
  )
) |>
  theme_latex(resize_direction = "down") |> 
  group_tt(j = list("Deflationary (2009-2021)" = 2:3, "Inflationary (2022-2023)" = 4:5))

save(profit.pre, profit.post, profit.pre.cov, profit.post.cov, profit.2023,
  profit.table,
  file = "stored_output/profit.prepost.RData"
)

## early placebo treatment estimation with 2004 -------------------------
set.seed(123)
early_data <- empdata |>
  mutate(
    treat = if_else(Bank == "BoE" & year >= 2004, 1, 0),
  )
class(early_data) <- "data.frame"
out.full.early <- bpCausal(
  data = early_data,
  Yname = "profit_ppgdp",
  Dname = "treat",
  Xname = c(),
  # Zname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Zname = c("traded", "Eurozone", "reappointable"),
  # Aname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Aname = c("traded", "Eurozone", "reappointable"),
  index = c("Bank", "year"),
  re = "both",
  ar1 = TRUE,
  r = 10,
  niter = 15000,
  burn = 5000,
  xlasso = 1,
  zlasso = 1,
  alasso = 1,
  flasso = 1,
  a1 = 0.001,
  a2 = 0.001,
  b1 = 0.001,
  b2 = 0.001,
  c1 = 0.001,
  c2 = 0.001,
  p1 = 0.001,
  p2 = 0.001
) # DM-LFM Warning: a bit time consuming
# see the plot
result.use.early <- effSummary(out.full.early)
placebo.plot <- ggplot(result.use.early$est.eff, aes(x = time + 2004)) +
  geom_ribbon(aes(ymin = counterfactual_ci_l, ymax = counterfactual_ci_u, fill = "95% CI"), alpha = 0.5) +
  geom_line(aes(y = observed, color = "Observed"), linewidth = 0.75, alpha = 0.5) +
  geom_line(aes(y = estimated_counterfactual, color = "Counterfactual"), linewidth = 0.75, alpha = 0.5) +
  geom_vline(xintercept = 2003, linetype = "dashed") +
  geom_vline(xintercept = 2008) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_color_manual(
    name = "",
    breaks = c("Observed", "Counterfactual"),
    values = c("Observed" = "#00BFC4", "Counterfactual" = "#F8766D")
  ) +
  labs(x = "Year", y = "BoE profit(% GDP)") +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(1999, 2023, 1), minor_breaks = NULL) + # ylim(-0.5, 0.5) +
  scale_fill_manual(
    name = c("", ""),
    breaks = c("95% CI"),
    values = c("lightgrey")
  )

profit.pre.early <- coefCustom(out.full.early, custom.period = 5 + (1:13))
profit.pre.early[1, ] <- effCustom(out.full.early, custom.period = 5 + (1:13))[1, ]
rownames(profit.pre.early) <- c("ATT", "traded", "Eurozone", "reappointable")

profit.post.early <- coefCustom(out.full.early, custom.period = 5 + (14:15))
profit.post.early[1, ] <- effCustom(out.full.early, custom.period = 5 + (14:15))[1, ]
rownames(profit.post.early) <- c("ATT", "traded", "Eurozone", "reappointable")

profit.pre.early <- list(results = profit.pre.early, obs = nrow(empdata), tr = 1, ct = length(unique(empdata$Country)) - 1)
profit.post.early <- list(results = profit.post.early, obs = nrow(empdata), tr = 1, ct = length(unique(empdata$Country)) - 1)
class(profit.pre.early) <- "DMLFM"
class(profit.post.early) <- "DMLFM"

save(out.full.early, result.use.early, profit.pre.early, profit.post.early, file = "stored_output/profit.early.RData")

cm <- c(
  "ATT" = "ATT",
  "traded" = "Publicly traded",
  "Eurozone" = "Euro Area",
  "reappointable" = "Reappointability"
)

profit.table.early <- modelsummary(
  list(
    "Early" = profit.pre.early,
    "Original" = profit.pre.cov,
    "Early" = profit.post.early,
    "Original" = profit.post.cov
  ),
  estimate = "{estimate}{stars}", output = "latex",
  statistic = "[{conf.low}, {conf.high}]",
  gof_map = gof_map, coef_map = cm, escape = F,
  notes = c(
    "+ p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001",
    "95\\% equal-tailed Credible Intervals in square brackets."
  )
) |>
  theme_latex(resize_direction = "down") |> 
  group_tt(j = list("Deflationary (2009-2021)" = 2:3, "Inflationary (2022-2023)" = 4:5))

save(profit.table.early, profit.pre.early, file = "stored_output/result.profit.early.RData")

# ir and balance-sheet dmlfm ---------------------------------------------
## ir dmlfm ----------------------------------------------------------------
ir_data <- empdata |>
  mutate(
    treat = if_else(Bank == "BoE" & year > 2008, 1, 0),
  ) |>
  filter(!is.na(ir))
class(ir_data) <- "data.frame"
set.seed(123)
out.ir <- bpCausal(
  data = ir_data,
  Yname = "ir",
  Dname = "treat",
  Xname = c(),
  Zname = c("traded", "Eurozone", "reappointable"),
  # Zname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Aname = c("traded", "Eurozone", "reappointable"),
  # Aname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  index = c("Bank", "year"),
  re = "both",
  ar1 = TRUE,
  r = 10,
  niter = 15000,
  burn = 5000,
  xlasso = 1,
  zlasso = 1,
  alasso = 1,
  flasso = 1,
  a1 = 0.001,
  a2 = 0.001,
  b1 = 0.001,
  b2 = 0.001,
  c1 = 0.001,
  c2 = 0.001,
  p1 = 0.001,
  p2 = 0.001
) # DM-LFM Warning: a bit time consuming
result.use.ir <- effSummary(out.ir) # extract effect information
ir.plot <- ggplot(result.use.ir$est.eff, aes(x = time + 2008)) +
  geom_ribbon(aes(ymin = counterfactual_ci_l, ymax = counterfactual_ci_u, fill = "95% CI"), alpha = 0.5) +
  geom_line(aes(y = observed, color = "Observed"), linewidth = 0.75, alpha = 0.5) +
  geom_line(aes(y = estimated_counterfactual, color = "Counterfactual"), linewidth = 0.75, alpha = 0.5) +
  geom_vline(xintercept = 2008) +
  geom_vline(xintercept = 2022, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_color_manual(
    name = "",
    breaks = c("Observed", "Counterfactual"),
    values = c("Observed" = "#00BFC4", "Counterfactual" = "#F8766D")
  ) +
  labs(x = "Year", y = "BoE policy interest rate (%)") +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(1999, 2023, 1), minor_breaks = NULL) +
  scale_fill_manual(
    name = c("", ""),
    breaks = c("95% CI"),
    values = c("lightgrey")
  )


### For previewing the trends only (not used in the paper)
ggplot(ir_data, aes(x = year, y = ir, color = Country)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Interest rate",
    subtitle = "%",
    x = "Year",
    y = "%Net gain/loss%"
  )

## balance sheet dmlfm -----------------------------------------------------
class(empdata) <- "data.frame"
set.seed(123)
out.bs <- bpCausal(
  data = empdata,
  # Yname = "profit_ppgdp",
  Yname = "liabilities_ppgdp",
  Dname = "treat",
  Xname = c(),
  Zname = c("traded", "Eurozone", "reappointable"),
  # Zname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  Aname = c("traded", "Eurozone", "reappointable"),
  # Aname = c('liabilities_ppgdp','ir','ir:liabilities_ppgdp'),
  index = c("Bank", "year"),
  re = "both",
  ar1 = TRUE,
  r = 10,
  niter = 15000,
  burn = 5000,
  xlasso = 1,
  zlasso = 1,
  alasso = 1,
  flasso = 1,
  a1 = 0.001,
  a2 = 0.001,
  b1 = 0.001,
  b2 = 0.001,
  c1 = 0.001,
  c2 = 0.001,
  p1 = 0.001,
  p2 = 0.001
) # DM-LFM Warning: a bit time consuming
result.use.bs <- effSummary(out.bs) # extract effect information
bs.plot <- ggplot(result.use.bs$est.eff, aes(x = time + 2008)) +
  geom_ribbon(aes(ymin = counterfactual_ci_l, ymax = counterfactual_ci_u, fill = "95% CI"), alpha = 0.5) +
  geom_line(aes(y = observed, color = "Observed"), linewidth = 0.75, alpha = 0.5) +
  geom_line(aes(y = estimated_counterfactual, color = "Counterfactual"), linewidth = 0.75, alpha = 0.5) +
  geom_vline(xintercept = 2008) +
  geom_vline(xintercept = 2022, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_color_manual(
    name = "",
    breaks = c("Observed", "Counterfactual"),
    values = c("Observed" = "#00BFC4", "Counterfactual" = "#F8766D")
  ) +
  labs(x = "Year", y = "BoE liabilities (% GDP)") +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(1999, 2023, 1), minor_breaks = NULL) +
  scale_fill_manual(
    name = c("", ""),
    breaks = c("95% CI"),
    values = c("lightgrey")
  )

empdata |>
  ggplot(aes(x = year, y = liabilities_ppgdp, color = Country)) +
  geom_line() +
  geom_point() +
  labs(
    title = "CB Liabilities",
    subtitle = "% GDP",
    x = "Year",
    y = "CB Liabilities (% GDP)"
  )

ir.pre <- coefCustom(out.ir, custom.period = 1:13)
ir.pre[1, ] <- effCustom(out.ir, custom.period = 1:13)[1, ]
rownames(ir.pre) <- c("ATT", "traded", "Eurozone", "reappointable")
ir.pre <- list(results = ir.pre, obs = nrow(ir_data), tr = 1, ct = length(unique(ir_data$Country)) - 1)
class(ir.pre) <- "DMLFM"

bs.pre <- coefCustom(out.bs, custom.period = 1:13)
bs.pre[1, ] <- effCustom(out.bs, custom.period = 1:13)[1, ]
rownames(bs.pre) <- c("ATT", "traded", "Eurozone", "reappointable")
bs.pre <- list(results = bs.pre, obs = nrow(empdata), tr = 1, ct = length(unique(empdata$Country)) - 1)
class(bs.pre) <- "DMLFM"

save(out.ir, ir.pre, file = "stored_output/supp-ir.RData")
save(out.bs, bs.pre, file = "stored_output/supp-bs.RData")

## supp table -------------------------------------------------------------
supp.table <- modelsummary(
  list(
    "Interest Rate" = ir.pre,
    "BoE Liabilities" = bs.pre
  ),
  estimate = "{estimate}{stars}", output = "latex",
  statistic = "[{conf.low}, {conf.high}]",
  gof_map = gof_map, coef_map = cm, escape = F,
  notes = c(
    "+ p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001",
    "95\\% equal-tailed Credible Intervals in square brackets."
  )
) |>
  theme_latex(resize_direction = "down")

save(supp.table, ir.plot, bs.plot, file = "stored_output/suppresults.RData")


# did/scm/sdid -------------------------------------------------------------
library(fixest)
library(synthdid)
sdid_setup <- empdata |>
  filter(!Country %in% c("Cyprus", "Japan", "Netherlands")) |>
  filter(year > 1999) |>
  mutate(treat = as.integer(treat)) |>
  as.data.frame() |>
  panel.matrices(unit = "Country", time = "year", outcome = "profit_ppgdp", treatment = "treat")

sdid_setup_pre <- empdata |>
  filter(!Country %in% c("Cyprus", "Japan", "Netherlands")) |>
  filter(year > 1999 & year < 2022) |>
  mutate(treat = as.integer(treat)) |>
  as.data.frame() |>
  panel.matrices(unit = "Country", time = "year", outcome = "profit_ppgdp", treatment = "treat")

sdid_setup_post <- empdata |>
  filter(!Country %in% c("Cyprus", "Japan", "Netherlands")) |>
  filter((year > 1999 & year < 2009) | year >= 2022) |>
  mutate(treat = as.integer(treat)) |>
  as.data.frame() |>
  panel.matrices(unit = "Country", time = "year", outcome = "profit_ppgdp", treatment = "treat")

did_data <- empdata |>
  filter(!Country %in% c("Cyprus", "Japan", "Netherlands")) |>
  filter(year > 1999) |>
  mutate(treat = as.integer(treat)) |>
  as.data.frame()

# synthdid_general for plotting
sdid <- synthdid_estimate(sdid_setup$Y, sdid_setup$N0, sdid_setup$T0)
did <- did_estimate(sdid_setup$Y, sdid_setup$N0, sdid_setup$T0)
scm <- sc_estimate(sdid_setup$Y, sdid_setup$N0, sdid_setup$T0)

# plot
source("custom_functions/solver.R")
synthdid_effect_custom <- function(estimate) {
  setup <- attr(estimate, "setup")
  weights <- attr(estimate, "weights")
  X.beta <- contract3(setup$X, weights$beta)
  N1 <- nrow(setup$Y) - setup$N0
  T1 <- ncol(setup$Y) - setup$T0

  tau.sc <- t(c(-weights$omega, rep(1 / N1, N1))) %*% (setup$Y - X.beta)
  # tau.curve = tau.sc[setup$T0 + (1:T1)] - c(tau.sc[1:setup$T0] %*% weights$lambda)
  tau.curve <- tau.sc[1:ncol(setup$Y)] - c(tau.sc[1:setup$T0] %*% weights$lambda)
  tau.curve
}
alt.plot.data <- data.frame(
  year = rep(2000:2023, 3),
  estimator = rep(c("SDiD", "DiD", "SCM"), each = 24),
  effect = c(synthdid_effect_custom(sdid), synthdid_effect_custom(did), synthdid_effect_custom(scm))
)
alt.plot.data <- rbind(
  alt.plot.data,
  data.frame(
    year = 1999:2023,
    estimator = rep("DM-LFM", each = 25),
    effect = result.use.full.cov$est.eff$estimated_ATT
  )
)

alt.plot <- ggplot(alt.plot.data, aes(x = year, y = effect, color = estimator)) +
  theme_bw() +
  geom_line() +
  geom_vline(xintercept = 2008) +
  geom_vline(xintercept = 2022, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_x_continuous(breaks = seq(1999, 2023, 1), minor_breaks = NULL) +
  labs(x = "Year", y = "Estimated effect on BoE profit (% GDP)", color = "Estimator") +
  theme(legend.position = "bottom")

## sdid -----------------------------------------------
sdid_pre <- synthdid_estimate(sdid_setup_pre$Y, sdid_setup_pre$N0, sdid_setup_pre$T0)
sdid_post <- synthdid_estimate(sdid_setup_post$Y, sdid_setup_post$N0, sdid_setup_post$T0)
set.seed(123)
attributes(sdid_pre)$se <- sqrt(vcov(sdid_pre, method = "placebo"))
attributes(sdid_post)$se <- sqrt(vcov(sdid_post, method = "placebo"))

alts <- list()
alts[["SDiD"]][["pre"]] <- data.frame(
  estimate = sdid_pre |> as.numeric(),
  p.value = 2 * pnorm(abs((sdid_pre |> as.numeric()) / attributes(sdid_pre)$se), lower.tail = FALSE)
)

alts[["SDiD"]][["post"]] <- data.frame(
  estimate = sdid_post |> as.numeric(),
  p.value = 2 * pnorm(abs((sdid_post |> as.numeric()) / attributes(sdid_post)$se), lower.tail = FALSE)
)

## did -----------------------------------------------

did_pre <- feols(profit_ppgdp ~ treat | Country + year, data = subset(did_data, year < 2022), cluster = "Country")
did_post <- feols(profit_ppgdp ~ treat | Country + year, data = subset(did_data, year >= 2022 | year < 2009), cluster = "Country")

alts[["DiD"]][["pre"]] <- get_estimates(did_pre)[, c(2, 9)]
alts[["DiD"]][["post"]] <- get_estimates(did_post)[, c(2, 9)]

## scm ----------------------------------------
scm_pre <- sc_estimate(sdid_setup_pre$Y, sdid_setup_pre$N0, sdid_setup_pre$T0)
scm_post <- sc_estimate(sdid_setup_post$Y, sdid_setup_post$N0, sdid_setup_post$T0)
source("custom_functions/RMSPE.R")

### permutation for pre-2022 ----------------------
scm_pre_rmspe <- scm_pre |> RMSPE_ratio()
for (i in 1:20) {
  sdid_setup_temp <- sdid_setup_pre
  sdid_setup_temp$Y[c(i, 21), ] <- sdid_setup_temp$Y[c(21, i), ]
  rownames(sdid_setup_temp$Y)[c(i, 21)] <- rownames(sdid_setup_temp$Y)[c(21, i)]
  rownames(sdid_setup_temp$W)[c(i, 21)] <- rownames(sdid_setup_temp$W)[c(21, i)]
  scm_temp <- sc_estimate(sdid_setup_temp$Y, sdid_setup_temp$N0, sdid_setup_temp$T0)
  scm_pre_rmspe <- append(scm_pre_rmspe, scm_temp |> RMSPE_ratio())
}

alts[["SCM"]][["pre"]] <- data.frame(
  estimate = scm_pre |> as.numeric(),
  p.value = sum(scm_pre_rmspe>=scm_pre_rmspe[1]) / length(scm_pre_rmspe)
)

### permutation for post-2022 ----------------------
scm_post_rmspe <- scm_post |> RMSPE_ratio()
for (i in 1:20) {
  sdid_setup_temp <- sdid_setup_post
  sdid_setup_temp$Y[c(i, 21), ] <- sdid_setup_temp$Y[c(21, i), ]
  rownames(sdid_setup_temp$Y)[c(i, 21)] <- rownames(sdid_setup_temp$Y)[c(21, i)]
  rownames(sdid_setup_temp$W)[c(i, 21)] <- rownames(sdid_setup_temp$W)[c(21, i)]
  scm_temp <- sc_estimate(sdid_setup_temp$Y, sdid_setup_temp$N0, sdid_setup_temp$T0)
  scm_post_rmspe <- append(scm_post_rmspe, scm_temp |> RMSPE_ratio())
}

alts[["SCM"]][["post"]] <- data.frame(
  estimate = scm_post |> as.numeric(),
  p.value = sum(scm_post_rmspe>=scm_post_rmspe[1]) / length(scm_post_rmspe)
)

## dm-lfm ------------------------------------------
alts[["DM-LFM"]][["pre"]] <- profit.pre.cov$results[1, c(1, 4)]
colnames(alts[["DM-LFM"]][["pre"]]) <- c("estimate", "p.value")
alts[["DM-LFM"]][["post"]] <- profit.post.cov$results[1, c(1, 4)]
colnames(alts[["DM-LFM"]][["post"]]) <- c("estimate", "p.value")

alts[["DM-LFM"]][["obs"]] <- profit.post.cov$obs
alts[["DM-LFM"]][["tr"]] <- profit.post.cov$tr
alts[["DM-LFM"]][["ct"]] <- profit.post.cov$ct

## numerical results ---------------------------------

for (i in 1:3) {
  alts[[i]][["obs"]] <- sdid_setup$Y |> length()
  alts[[i]][["tr"]] <- 1
  alts[[i]][["ct"]] <- sdid_setup$N0
}

for (i in 1:4) class(alts[[i]]) <- "alt"

### table -----------------------------------------
tidy.alt <- function(x, ...) {
  ret <- data.frame(
    term      = c("pre", "post"),
    estimate  = c(x$pre$estimate, x$post$estimate),
    p.value   = c(x$pre$p.value, x$post$p.value)
  )
  ret
}

glance.alt <- function(x, ...) {
  ret <- data.frame(
    nobs = x$obs,
    tr = x$tr,
    ct = x$ct
  )
  ret
}

gof_map <- list(
  list("raw" = "nobs", "clean" = "Observations", "fmt" = 0),
  list("raw" = "tr", "clean" = "Treated Units", "fmt" = 0),
  list("raw" = "ct", "clean" = "Control Units", "fmt" = 0)
)
cm <- c(
  "pre" = "$ATT^\\text{def}$",
  "post" = "$ATT^\\text{inf}$"
)

alt.table <- modelsummary(alts,
  output = "latex", # comment this out for testing
  statistic = "p.value",
  stars = T,
  gof_map = gof_map, coef_map = cm, escape = F,
  notes = "$p$ values are reported in parentheses."
)

save(alt.table, alt.plot, file = "stored_output/alt.RData")
