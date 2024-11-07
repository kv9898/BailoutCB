need <- c("bayestestR") # list packages needed
have <- need %in% rownames(installed.packages()) # checks packages you have
if(any(!have)) install.packages(need[!have]) # install missing packages
#devtools::install_github("synth-inference/synthdid") #manually install sdid if the codes above fail to do so automatically
invisible(lapply(need, library, character.only=T)) # load needed packages

effSummary <- function(x, 
                       usr.id = NULL,         ## individual effect, if left blank, all treated units will be used
                       burn = 0, 
                       cumu = FALSE,          ## whether to calculate cumulative effect
                       rela.period = TRUE,    ## aggregate by time relative to treatment
                       ci.level = 0.95,       ## set CI level (not limited to 0.95)
                       ci.type  = TRUE) {     ## TRUE for equal-tailed CI, or FALSE for "HDI"
  
  niter <- dim(x$sigma2)[2] 
  
  if (cumu) {
    rela.period <- TRUE
  }
  
  id.tr <- x$raw.id.tr
  time.tr <- x$time.tr
  rela.time.tr <- x$rela.time.tr  
  
  
  ## id indicator
  id.pos <- NULL
  unique.tr <- c(unique(id.tr))
  if (is.null(usr.id)) {
    id.pos <- 1:length(c(id.tr))  
  } else {
    if (sum(usr.id %in% unique.tr) != length(usr.id)) {
      stop("Some specified ids are not in treated group, please check input.\n")
    }
    id.pos <- which(c(id.tr) %in% usr.id)
  }
  
  yo_t <- x$yo_t
  yo_t <- yo_t[id.pos]
  
  time.tr <- time.tr[id.pos]
  rela.time.tr <- rela.time.tr[id.pos]
  
  yct_i <- x$yct
  yct_i <- matrix(c(yct_i[id.pos, (burn + 1):niter]), length(id.pos), niter - burn)
  
  ## mean observed and counterfactual
  count.tr <- NULL ## num of observations at each period
  
  if (rela.period) { ## relative to treatment occurence
    m_yo <- tapply(yo_t, rela.time.tr, mean)
    m_yct <- sapply(1:(niter - burn), function(i){tapply(yct_i[, i], rela.time.tr, mean)})
    
    count.tr <- as.numeric(table(rela.time.tr))
  } else { ## real time period
    m_yo <- tapply(yo_t, time.tr, mean)
    m_yct <- sapply(1:(niter - burn), function(i){tapply(yct_i[, i], time.tr, mean)})
    
    count.tr <- as.numeric(table(rela.time.tr))
  }
  
  
  ## outcomes -------------------
  
  m_yct_mean <- apply(m_yct, 1, mean)
  m_yct_ci_l <- if(ci.type) apply(m_yct, 1, quantile, (1-ci.level)/2) else apply(m_yct, 1, function(x, ci) hdi(x,ci)$CI_low, ci=ci.level)
  m_yct_ci_u <- if(ci.type) apply(m_yct, 1, quantile, 1-((1-ci.level)/2)) else apply(m_yct, 1, function(x, ci) hdi(x,ci)$CI_high, ci=ci.level)
  
  
  ## effect ---------------------
  
  eff_i <- matrix(rep(c(m_yo), niter - burn), length(c(m_yo)), niter - burn) - m_yct
  
  eff_mean <- apply(eff_i, 1, mean)
  eff_ci_l <- if(ci.type) apply(eff_i, 1, quantile, (1-ci.level)/2) else apply(eff_i, 1, function(x, ci) hdi(x,ci)$CI_low, ci=ci.level)
  eff_ci_u <- if(ci.type) apply(eff_i, 1, quantile, 1-((1-ci.level)/2)) else apply(eff_i, 1, function(x, ci) hdi(x,ci)$CI_high, ci=ci.level)
  
  
  data <- cbind.data.frame(m_yo, m_yct_mean, m_yct_ci_u, m_yct_ci_l, eff_mean, eff_ci_l, eff_ci_u)
  names(data) <- c("observed", "estimated_counterfactual", 
                   "counterfactual_ci_l", "counterfactual_ci_u",
                   "estimated_ATT", "estimated_ATT_ci_l", "estimated_ATT_ci_u")
  if(rela.period) {
    data$time <- sort(unique(rela.time.tr))
    data$count <- count.tr
  } else {
    data$time <- sort(unique(time.tr))
  }
  
  est.eff <- data
  
  
  ## cumulative effects ---------
  est.cumu <- NULL
  if (cumu) {
    relatime <- sort(unique(rela.time.tr))
    
    st.pos <- which(relatime == 1) ## start point 
    
    eff_sub_i <- matrix(c(eff_i[st.pos:length(relatime), ]), length(relatime) - st.pos +1)
    
    eff_cumu_i <- matrix(NA, length(relatime) - st.pos +1, niter - burn)
    
    count.tr.sub <- count.tr[st.pos:length(relatime)]
    
    eff_cumu_i[1, ] <- eff_sub_i[1, ]
    if (length(relatime) - st.pos >= 2) {
      for (j in 2:(length(relatime) - st.pos +1)) {
        eff_cumu_i[j, ] <- sapply(1:(niter - burn), function(i) {sum(eff_sub_i[1:j, i] * count.tr.sub[1:j])/sum(count.tr.sub[1:j])} ) * j
      }
    }
    
    eff_cumu_mean <- apply(eff_cumu_i, 1, mean)
    eff_cumu_ci_l <- if(ci.type) apply(eff_cumu_i, 1, quantile, (1-ci.level)/2) else apply(eff_cumu_i, 1, function(x, ci) hdi(x,ci)$CI_low, ci=ci.level)
    eff_cumu_ci_u <- if(ci.type) apply(eff_cumu_i, 1, quantile, 1-((1-ci.level)/2)) else apply(eff_cumu_i, 1, function(x, ci) hdi(x,ci)$CI_high, ci=ci.level)
    
    data <- cbind.data.frame(eff_cumu_mean, eff_cumu_ci_l, eff_cumu_ci_u)
    names(data) <- c("mean", "ci_l", "ci_u")
    
    data$count <- count.tr[st.pos:length(relatime)]
    data$time <- relatime[st.pos:length(relatime)]
    
    est.cumu <- data
  }
  
  
  
  ## average effects ------------
  
  t.post <- which(rela.time.tr > 0)
  
  eff_avg_i <- sapply(1:(niter - burn), function(i) {mean(yo_t[t.post] - yct_i[t.post, i])})
  
  eff_avg_mean <- mean(eff_avg_i)
  eff_avg_ci_l <- if(ci.type) quantile(eff_avg_i, (1-ci.level)/2) else hdi(eff_avg_i, ci=ci.level)$CI_low
  eff_avg_ci_u <- if(ci.type) quantile(eff_avg_i, 1-((1-ci.level)/2)) else hdi(eff_avg_i, ci=ci.level)$CI_high
  
  #p-value calculation
  eff_avg_p <- sum(eff_avg_i<=0)/length(eff_avg_i)
  eff_avg_p <- ifelse(eff_avg_p<0.5, eff_avg_p*2, (1-eff_avg_p)*2)
  
  
  est.avg <- cbind(eff_avg_mean, eff_avg_ci_l, eff_avg_ci_u)
  colnames(est.avg) <- c("mean", "ci_l", "ci_u")
  
  
  out <- list(est.eff = est.eff, 
              est.avg = est.avg,
              est.avg_p = eff_avg_p,
              eff_avg_i = eff_avg_i)
  
  if (!is.null(est.cumu)) {
    out <- c(out, list(est.cumu = est.cumu))
  }
  
  return(out)
  
}

plot.DMLFM <- function(x, 
                       usr.id = NULL,         ## individual effect, if left blank, all treated units will be used
                       burn = 0, 
                       cumu = FALSE,          ## whether to calculate cumulative effect
                       rela.period = TRUE,    ## aggregate by time relative to treatment
                       ci.level = 0.95,       ## set CI level (not limited to 0.95)
                       ci.type  = TRUE) {     ## TRUE for equal-tailed CI, or FALSE for "HDI"
  niter <- dim(x$sigma2)[2] 
  
  if (cumu) {
    rela.period <- TRUE
  }
  
  id.tr <- x$raw.id.tr
  time.tr <- x$time.tr
  rela.time.tr <- x$rela.time.tr  
  
  
  ## id indicator
  id.pos <- NULL
  unique.tr <- c(unique(id.tr))
  if (is.null(usr.id)) {
    id.pos <- 1:length(c(id.tr))  
  } else {
    if (sum(usr.id %in% unique.tr) != length(usr.id)) {
      stop("Some specified ids are not in treated group, please check input.\n")
    }
    id.pos <- which(c(id.tr) %in% usr.id)
  }
  
  yo_t <- x$yo_t
  yo_t <- yo_t[id.pos]
  
  time.tr <- time.tr[id.pos]
  rela.time.tr <- rela.time.tr[id.pos]
  
  yct_i <- x$yct
  yct_i <- matrix(c(yct_i[id.pos, (burn + 1):niter]), length(id.pos), niter - burn)
  
  ## mean observed and counterfactual
  count.tr <- NULL ## num of observations at each period
  
  if (rela.period) { ## relative to treatment occurence
    m_yo <- tapply(yo_t, rela.time.tr, mean)
    m_yct <- sapply(1:(niter - burn), function(i){tapply(yct_i[, i], rela.time.tr, mean)})
    
    count.tr <- as.numeric(table(rela.time.tr))
  } else { ## real time period
    m_yo <- tapply(yo_t, time.tr, mean)
    m_yct <- sapply(1:(niter - burn), function(i){tapply(yct_i[, i], time.tr, mean)})
    
    count.tr <- as.numeric(table(rela.time.tr))
  }
  
  
  ## outcomes -------------------
  
  m_yct_mean <- apply(m_yct, 1, mean)
  m_yct_ci_l <- if(ci.type) apply(m_yct, 1, quantile, (1-ci.level)/2) else apply(m_yct, 1, function(x, ci) hdi(x,ci)$CI_low, ci=ci.level)
  m_yct_ci_u <- if(ci.type) apply(m_yct, 1, quantile, 1-((1-ci.level)/2)) else apply(m_yct, 1, function(x, ci) hdi(x,ci)$CI_high, ci=ci.level)
  
  
  ## effect ---------------------
  
  eff_i <- matrix(rep(c(m_yo), niter - burn), length(c(m_yo)), niter - burn) - m_yct
  
  eff_mean <- apply(eff_i, 1, mean)
  eff_ci_l <- if(ci.type) apply(eff_i, 1, quantile, (1-ci.level)/2) else apply(eff_i, 1, function(x, ci) hdi(x,ci)$CI_low, ci=ci.level)
  eff_ci_u <- if(ci.type) apply(eff_i, 1, quantile, 1-((1-ci.level)/2)) else apply(eff_i, 1, function(x, ci) hdi(x,ci)$CI_high, ci=ci.level)
  
  
  data <- cbind.data.frame(m_yo, m_yct_mean, m_yct_ci_u, m_yct_ci_l, eff_mean, eff_ci_l, eff_ci_u)
  names(data) <- c("observed", "estimated_counterfactual", 
                   "counterfactual_ci_l", "counterfactual_ci_u",
                   "estimated_ATT", "estimated_ATT_ci_l", "estimated_ATT_ci_u")
  est.eff <- data
  
  
  ## cumulative effects ---------
  est.cumu <- NULL
  if (cumu) {
    relatime <- sort(unique(rela.time.tr))
    
    st.pos <- which(relatime == 1) ## start point 
    
    eff_sub_i <- matrix(c(eff_i[st.pos:length(relatime), ]), length(relatime) - st.pos +1)
    
    eff_cumu_i <- matrix(NA, length(relatime) - st.pos +1, niter - burn)
    
    count.tr.sub <- count.tr[st.pos:length(relatime)]
    
    eff_cumu_i[1, ] <- eff_sub_i[1, ]
    if (length(relatime) - st.pos >= 2) {
      for (j in 2:(length(relatime) - st.pos +1)) {
        eff_cumu_i[j, ] <- sapply(1:(niter - burn), function(i) {sum(eff_sub_i[1:j, i] * count.tr.sub[1:j])/sum(count.tr.sub[1:j])} ) * j
      }
    }
    
    eff_cumu_mean <- apply(eff_cumu_i, 1, mean)
    eff_cumu_ci_l <- if(ci.type) apply(eff_cumu_i, 1, quantile, (1-ci.level)/2) else apply(eff_cumu_i, 1, function(x, ci) hdi(x,ci)$CI_low, ci=ci.level)
    eff_cumu_ci_u <- if(ci.type) apply(eff_cumu_i, 1, quantile, 1-((1-ci.level)/2)) else apply(eff_cumu_i, 1, function(x, ci) hdi(x,ci)$CI_high, ci=ci.level)
    
    data <- cbind.data.frame(eff_cumu_mean, eff_cumu_ci_l, eff_cumu_ci_u)
    names(data) <- c("mean", "ci_l", "ci_u")
    
    data$count <- count.tr[st.pos:length(relatime)]
    data$time <- relatime[st.pos:length(relatime)]
    
    est.cumu <- data
  }
  
  
  
  ## average effects ------------
  
  t.post <- which(rela.time.tr > 0)
  
  eff_avg_i <- sapply(1:(niter - burn), function(i) {mean(yo_t[t.post] - yct_i[t.post, i])})
  
  eff_avg_mean <- mean(eff_avg_i)
  eff_avg_ci_l <- if(ci.type) quantile(eff_avg_i, (1-ci.level)/2) else hdi(eff_avg_i, ci=ci.level)$CI_low
  eff_avg_ci_u <- if(ci.type) quantile(eff_avg_i, 1-((1-ci.level)/2)) else hdi(eff_avg_i, ci=ci.level)$CI_high
  
  p_test <- function(value){
    eff_avg_p <- sapply(value, function(value) {
      p_value <- sum(eff_avg_i <= value) / length(eff_avg_i)
      p_value <- ifelse(p_value < 0.5, p_value * 2, (1 - p_value) * 2)
      return(p_value)
    })
    return(eff_avg_p)
  }
  plot_data <- data.frame(x=eff_avg_i, y=p_test(eff_avg_i))
  ggplot(plot_data, aes(x=x, y=y))+geom_line()+geom_vline(xintercept=eff_avg_mean, colour='#00BFC4')+
    geom_vline(xintercept=eff_avg_ci_l,colour='blue',alpha=0.5)+
    geom_vline(xintercept=eff_avg_ci_u,colour='blue',alpha=0.5)+
    geom_hline(yintercept=p_test(0), colour='#F8766D')+labs(x='ATT',y='p-value')+
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+theme_bw()
  
  }

effCustom <- function(x, 
    custom.period = NULL,      ## custom treatment period for ATT calculation. If left blank, all treated periods will be used
    usr.id = NULL,         ## individual effect, if left blank, all treated units will be used
    burn = 0, 
    rela.period = TRUE,    ## aggregate by time relative to treatment
    ci.level = 0.95,       ## set CI level (not limited to 0.95)
    ci.type  = TRUE) {     ## TRUE for equal-tailed CI, or FALSE for "HDI"
  
    niter <- dim(x$sigma2)[2] 
    
    
    id.tr <- x$raw.id.tr
    time.tr <- x$time.tr
    rela.time.tr <- x$rela.time.tr  
    
    
    ## id indicator
    id.pos <- NULL
    unique.tr <- c(unique(id.tr))
    if (is.null(usr.id)) {
      id.pos <- 1:length(c(id.tr))  
    } else {
      if (sum(usr.id %in% unique.tr) != length(usr.id)) {
        stop("Some specified ids are not in treated group, please check input.\n")
      }
      id.pos <- which(c(id.tr) %in% usr.id)
    }
    
    yo_t <- x$yo_t
    yo_t <- yo_t[id.pos]
    
    time.tr <- time.tr[id.pos]
    rela.time.tr <- rela.time.tr[id.pos]
    
    yct_i <- x$yct
    yct_i <- matrix(c(yct_i[id.pos, (burn + 1):niter]), length(id.pos), niter - burn)
  
    if(is.null(custom.period)) t.post <- which(rela.time.tr > 0) else t.post <- which(rela.time.tr %in% custom.period)
    
  eff_avg_i <- sapply(1:(niter - burn), function(i) {mean(yo_t[t.post] - yct_i[t.post, i])})
  
  eff_avg_mean <- mean(eff_avg_i)
  eff_avg_ci_l <- if(ci.type) quantile(eff_avg_i, (1-ci.level)/2) else hdi(eff_avg_i, ci=ci.level)$CI_low
  eff_avg_ci_u <- if(ci.type) quantile(eff_avg_i, 1-((1-ci.level)/2)) else hdi(eff_avg_i, ci=ci.level)$CI_high
  
  #p-value calculation
  eff_avg_p <- sum(eff_avg_i<=0)/length(eff_avg_i)
  eff_avg_p <- ifelse(eff_avg_p<0.5, eff_avg_p*2, (1-eff_avg_p)*2)
  
  est.avg <- cbind(eff_avg_mean, eff_avg_ci_l, eff_avg_ci_u, eff_avg_p)
    colnames(est.avg) <- c("mean", "ci_l", "ci_u", "p")
    
    
    return(est.avg)
  
    }

    coefCustom <- function(x,                 ## estimation results
      custom.period = NULL,      ## custom treatment period for ATT calculation. If left blank, all treated periods will be used
      ci.level = 0.95,       ## set CI level (not limited to 0.95)
      ci.type  = TRUE,       ## TRUE for equal-tailed CI, or FALSE for "HDI"
      burn = 0) {        ## burn-in length
    
    niter <- dim(x$sigma2)[2]
    if(is.null(custom.period)) t.post <- which(x$rela.time.tr > 0) else t.post <- which(x$rela.time.tr %in% custom.period)
    
    ## 3. time-varying coefficient 
    
    xi_i <- x$xi 
    wxi_i <- x$wxi
    
    if (!is.null(xi_i)) {
    TT <- dim(xi_i)[1]
    p <- dim(xi_i)[2]
    
    for (i in (burn + 1):niter) {
    xi_i[,, i] <- xi_i[,,i] * matrix(rep(wxi_i[, i], each = TT), TT, p)
    }
    
    est_xi_ci_l <- est_xi_ci_u <- est_xi_mean <- matrix(NA, TT, p)
      
    output <- data.frame(mean = numeric(), ci_l = numeric(), ci_u = numeric(), p=numeric())
      
    for ( i in 1:p){
      sub_xi <- matrix(xi_i[t.post,i,], length(t.post), niter)
      est_xi_i <- colMeans(sub_xi)
      output[i,1] <- mean(est_xi_i)
      output[i,2] <- if(ci.type) quantile(est_xi_i, (1-ci.level)/2) else hdi(est_xi_i, ci=ci.level)$CI_low
      output[i,3] <- if(ci.type) quantile(est_xi_i, 1-((1-ci.level)/2)) else hdi(est_xi_i, ci=ci.level)$CI_high
      p <- sum(est_xi_i<=0)/length(est_xi_i)
      output[i,4] <-  ifelse(p<0.5, p*2, (1-p)*2)
    }
    
    }
    return(output)
    } 
    