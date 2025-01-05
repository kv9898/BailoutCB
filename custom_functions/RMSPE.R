source("custom_functions/solver.R")
RMSPE_ratio <- function (estimate) 
{
  setup = attr(estimate, "setup")
  weights = attr(estimate, "weights")
  X.beta = contract3(setup$X, weights$beta)
  N1 = nrow(setup$Y) - setup$N0
  T1 = ncol(setup$Y) - setup$T0
  tau.sc = t(c(-weights$omega, rep(1/N1, N1))) %*% (setup$Y - 
                                                      X.beta)
  tau.pre.curve = tau.sc[1:setup$T0] - c(tau.sc[1:setup$T0] %*% 
                                           weights$lambda)
  tau.post.curve = tau.sc[setup$T0+(1:T1)] - c(tau.sc[1:setup$T0] %*% 
                                                 weights$lambda)
  RMSPE.pre <- sqrt(sum(tau.pre.curve^2)/length(tau.pre.curve))
  RMSPE.post <- sqrt(sum(tau.post.curve^2)/length(tau.post.curve))
  RMSPE.post/RMSPE.pre
}