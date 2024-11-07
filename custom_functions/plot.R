source("custom_functions/solver.r")
synthdid_plot <- function (estimates, treated.name = "treated", control.name = "synthetic control", 
          spaghetti.units = c(), spaghetti.matrices = NULL, facet = NULL, 
          facet.vertical = TRUE, lambda.comparable = !is.null(facet), 
          overlay = 0, lambda.plot.scale = 3, trajectory.linetype = 1, 
          effect.curvature = 0.3, line.width = 0.5, guide.linetype = 2, 
          point.size = 1, trajectory.alpha = 0.5, diagram.alpha = 0.95, 
          effect.alpha = 0.95, onset.alpha = 0.3, ci.alpha = 0.3, 
          spaghetti.line.width = 0.2, spaghetti.label.size = 2, spaghetti.line.alpha = 0.3, 
          spaghetti.label.alpha = 0.5, se.method = "jackknife", se.manual=NULL, ci.manual=NULL, alpha.multiplier = NULL) 
{
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .ignore <- tryCatch(attachNamespace("ggplot2"), error = function(e) e)
  }
  else {
    stop("Plotting requires the package `ggplot2`. Install it to use this function.")
  }
  if (class(estimates) == "synthdid_estimate") {
    estimates = list(estimates)
  }
  if (is.null(names(estimates))) {
    names(estimates) = sprintf("estimate %d", 1:length(estimates))
  }
  if (is.null(alpha.multiplier)) {
    alpha.multiplier = rep(1, length(estimates))
  }
  if (!is.null(spaghetti.matrices) && length(spaghetti.matrices) != 
      length(estimates)) {
    stop("spaghetti.matrices must be the same length as estimates")
  }
  multiple.frames = length(overlay) > 1
  treated = 1
  control = 2
  groups = factor(c(control, treated), labels = c(control.name, 
                                                  treated.name))
  estimate.factors = factor(1:(length(estimates) + 1), labels = c(treated.name, 
                                                                  names(estimates)))
  facet_factors = if (is.null(facet)) {
    factor(1:length(estimates), labels = names(estimates))
  }
  else {
    factor(facet, levels = 1:length(unique(facet)), labels = unique(facet))
  }
  grid = expand.grid(estimate = 1:length(estimates), overlay = 1:length(overlay))
  plot.descriptions = lapply(1:nrow(grid), function(row) {
    est = estimates[[grid$estimate[row]]]
    over = overlay[grid$overlay[row]]
    se = if (se.method == "none") {
      NA
    }
    else {
      if (!is.null(se.manual)){
        se.manual[[row]]
      } else{
        sqrt(vcov(est, method = se.method))
      }
      
    }
    setup = attr(est, "setup")
    weights = attr(est, "weights")
    Y = setup$Y - contract3(setup$X, weights$beta)
    N0 = setup$N0
    N1 = nrow(Y) - N0
    T0 = setup$T0
    T1 = ncol(Y) - T0
    lambda.synth = c(weights$lambda, rep(0, T1))
    lambda.target = c(rep(0, T0), rep(1/T1, T1))
    omega.synth = c(weights$omega, rep(0, N1))
    omega.target = c(rep(0, N0), rep(1/N1, N1))
    if (!is.null(attr(est, "overlay"))) {
      over = attr(est, "overlay")
    }
    is.sc = all(weights$lambda == 0) || over == 1
    intercept.offset = over * c((omega.target - omega.synth) %*% 
                                  Y %*% lambda.synth)
    obs.trajectory = as.numeric(omega.target %*% Y)
    syn.trajectory = as.numeric(omega.synth %*% Y) + intercept.offset
    spaghetti.trajectories = Y[rownames(Y) %in% spaghetti.units, 
                               , drop = FALSE]
    if (!is.null(spaghetti.matrices)) {
      more.spaghetti.trajectories = spaghetti.matrices[[grid$estimate[row]]]
      if (ncol(more.spaghetti.trajectories) != ncol(Y)) {
        stop("The elements of spaghetti.matrices must be matrices with the same number of columns as Y")
      }
      if (is.null(rownames(more.spaghetti.trajectories))) {
        stop("The elements of the list spaghetti.matrices must have named rows")
      }
      spaghetti.trajectories = rbind(spaghetti.trajectories, 
                                     more.spaghetti.trajectories)
    }
    treated.post = omega.target %*% Y %*% lambda.target
    treated.pre = omega.target %*% Y %*% lambda.synth
    control.post = omega.synth %*% Y %*% lambda.target + 
      intercept.offset
    control.pre = omega.synth %*% Y %*% lambda.synth + intercept.offset
    sdid.post = as.numeric(control.post + treated.pre - 
                             control.pre)
    time = as.numeric(timesteps(Y))
    if (length(time) == 0 || !all(is.finite(time))) {
      time = 1:(T0 + T1)
    }
    pre.time = lambda.synth %*% time
    post.time = lambda.target %*% time
    lines = data.frame(x = rep(time, 2), y = c(obs.trajectory, 
                                               syn.trajectory), color = rep(groups[c(treated, control)], 
                                                                            each = length(time)))
    points = data.frame(x = c(post.time, post.time), y = c(treated.post, 
                                                           sdid.post), color = groups[c(treated, control)])
    did.points = data.frame(x = c(pre.time, pre.time, post.time, 
                                  post.time), y = c(treated.pre, control.pre, control.post, 
                                                    treated.post), color = groups[c(treated, control, 
                                                                                    control, treated)])
    did.segments = data.frame(x = c(pre.time, pre.time), 
                              xend = c(post.time, post.time), y = c(control.pre, 
                                                                    treated.pre), yend = c(control.post, treated.post), 
                              color = groups[c(control, treated)])
    hallucinated.segments = data.frame(x = pre.time, xend = post.time, 
                                       y = treated.pre, yend = sdid.post)
    guide.segments = data.frame(x = c(pre.time, post.time), 
                                xend = c(pre.time, post.time), y = c(control.pre, 
                                                                     control.post), yend = c(treated.pre, sdid.post))
    arrows = data.frame(x = post.time, xend = post.time, 
                        y = sdid.post, yend = treated.post, xscale = max(time) - 
                          post.time, color = groups[control])
    ub.arrows = data.frame(x = post.time, xend = post.time, 
                           y = ifelse(!is.null(ci.manual),treated.post-ci.manual[[row*2-1]],sdid.post + 1.96 * se), yend = treated.post, 
                           xscale = max(time) - post.time, color = groups[control])
    lb.arrows = data.frame(x = post.time, xend = post.time, 
                           y = ifelse(!is.null(ci.manual),treated.post-ci.manual[[row*2]],sdid.post - 1.96 * se), yend = treated.post, 
                           xscale = max(time) - post.time, color = groups[control])
    spaghetti.lines = data.frame(x = rep(time, nrow(spaghetti.trajectories)), 
                                 y = as.vector(t(spaghetti.trajectories)), unit = rep(rownames(spaghetti.trajectories), 
                                                                                      each = length(time)))
    spaghetti.labels = data.frame(x = rep(time[1], nrow(spaghetti.trajectories)), 
                                  y = as.vector(spaghetti.trajectories[, 1]), unit = rownames(spaghetti.trajectories))
    T0s = attr(est, "T0s")
    if (!is.null(T0s)) {
      vlines = data.frame(xintercept = time[T0s])
    }
    else {
      vlines = data.frame(xintercept = time[T0])
    }
    if (lambda.comparable) {
      height = (max(c(obs.trajectory)) - min(c(obs.trajectory)))/lambda.plot.scale
      bottom = min(c(obs.trajectory)) - height
      ribbons = data.frame(x = time[1:T0], ymin = rep(bottom, 
                                                      T0), ymax = bottom + height * lambda.synth[1:T0], 
                           color = groups[control])
    }
    else {
      height = (max(c(obs.trajectory, syn.trajectory)) - 
                  min(c(obs.trajectory, syn.trajectory)))/lambda.plot.scale
      bottom = min(c(obs.trajectory, syn.trajectory)) - 
        height
      ribbons = data.frame(x = time[1:T0], ymin = rep(bottom, 
                                                      T0), ymax = bottom + height * lambda.synth[1:T0]/max(lambda.synth), 
                           color = groups[control])
    }
    elements = list(lines = lines, points = points, did.segments = did.segments, 
                    did.points = did.points, hallucinated.segments = hallucinated.segments, 
                    guide.segments = guide.segments, arrows = arrows, 
                    lb.arrows = lb.arrows, ub.arrows = ub.arrows, spaghetti.lines = spaghetti.lines, 
                    spaghetti.labels = spaghetti.labels, vlines = vlines, 
                    ribbons = ribbons)
    lapply(elements, function(x) {
      if (nrow(x) > 0) {
        x$frame = over
        x$is.sc = is.sc
        x$estimate = estimate.factors[grid$estimate[row] + 
                                        1]
      }
      x
    })
  })
  one.per.facet = length(unique(facet_factors)) == length(facet_factors)
  concatenate.field = function(field) {
    do.call(rbind, lapply(plot.descriptions, function(desc) {
      element = desc[[field]]
      estimate.factor = element$estimate[1]
      element$facet = facet_factors[as.integer(estimate.factor) - 
                                      1]
      element$show = alpha.multiplier[as.integer(element$estimate) - 
                                        1]
      element$show[element$color == groups[treated]] = 1
      if (!one.per.facet && "color" %in% colnames(element)) {
        color = element$estimate
        color[element$color == groups[treated]] = estimate.factors[1]
        element$color = color
      }
      element
    }))
  }
  conc = lapply(names(plot.descriptions[[1]]), concatenate.field)
  names(conc) = names(plot.descriptions[[1]])
  no.sc = function(x) {
    x[!x$is.sc, ]
  }
  with.frame = function(geom, base.aes, data, ...) {
    new.aes = if (multiple.frames) {
      modifyList(base.aes, aes(frame = frame))
    }
    else {
      base.aes
    }
    do.call(geom, c(list(new.aes, data = data), list(...)))
  }
  p = ggplot() + with.frame(geom_line, aes(x = x, y = y, color = color, 
                                           alpha = trajectory.alpha * show), data = conc$lines, 
                            linetype = trajectory.linetype, linewidth = line.width) + 
    with.frame(geom_point, aes(x = x, y = y, color = color, 
                               alpha = diagram.alpha * show), data = conc$points, 
               shape = 21, size = point.size) + with.frame(geom_point, 
                                                           aes(x = x, y = y, color = color, alpha = diagram.alpha * 
                                                                 show), data = no.sc(conc$did.points), size = point.size) + 
    with.frame(geom_segment, aes(x = x, xend = xend, y = y, 
                                 yend = yend, color = color, alpha = diagram.alpha * 
                                   show), data = no.sc(conc$did.segments), size = line.width) + 
    with.frame(geom_segment, aes(x = x, xend = xend, y = y, 
                                 yend = yend, group = estimate, alpha = 0.6 * diagram.alpha * 
                                   show), data = no.sc(conc$hallucinated.segments), 
               linetype = guide.linetype, size = line.width, color = "black") + 
    with.frame(geom_segment, aes(x = x, xend = xend, y = y, 
                                 yend = yend, group = estimate, alpha = 0.5 * diagram.alpha * 
                                   show), data = no.sc(conc$guide.segments), size = line.width, 
               linetype = guide.linetype, color = "black") + geom_vline(aes(xintercept = xintercept, 
                                                                            alpha = onset.alpha * show), data = conc$vlines, size = line.width, 
                                                                        color = "black") + geom_ribbon(aes(x = x, ymin = ymin, 
                                                                                                           ymax = ymax, group = color, fill = color, alpha = 0.5 * 
                                                                                                             diagram.alpha * show), data = no.sc(conc$ribbons), 
                                                                                                       color = "black", size = line.width, show.legend = FALSE) + 
    geom_curve(aes(x = x, xend = xend, y = y, yend = yend, 
                   alpha = effect.alpha * show), data = conc$arrows, 
               curvature = effect.curvature, color = "black", size = line.width, 
               arrow = arrow(length = unit(0.2, "cm"))) + geom_curve(aes(x = x, 
                                                                         xend = xend, y = y, yend = yend, alpha = ci.alpha * 
                                                                           show), data = conc$ub.arrows, na.rm = TRUE, curvature = effect.curvature, 
                                                                     color = "black", size = line.width, arrow = arrow(length = unit(0.2, 
                                                                                                                                     "cm"))) + geom_curve(aes(x = x, xend = xend, y = y, 
                                                                                                                                                              yend = yend, alpha = ci.alpha * show), data = conc$lb.arrows, 
                                                                                                                                                          na.rm = TRUE, curvature = effect.curvature, color = "black", 
                                                                                                                                                          size = line.width, arrow = arrow(length = unit(0.2, 
                                                                                                                                                                                                         "cm")))
  if (nrow(conc$spaghetti.labels) > 0) {
    p = p + geom_text(aes(x = x, y = y, label = unit, alpha = spaghetti.label.alpha * 
                            show), data = conc$spaghetti.labels, color = "black", 
                      size = spaghetti.label.size) + geom_line(aes(x = x, 
                                                                   y = y, group = unit, alpha = spaghetti.line.alpha * 
                                                                     show), data = conc$spaghetti.lines, color = "black", 
                                                               size = spaghetti.line.width)
  }
  if (!all(conc$lines$facet == conc$lines$facet[1])) {
    if (facet.vertical) {
      p = p + facet_grid(facet ~ ., scales = "free_y")
    }
    else {
      p = p + facet_grid(. ~ facet)
    }
  }
  if (is.null(facet)) {
    p = p + guides(linetype = "none")
  }
  p = tryCatch({
    as.Date(colnames(attr(estimates[[1]], "setup")$Y))
    p + scale_x_continuous(labels = function(time) {
      as.Date(time, origin = "1970-01-01")
    })
  }, error = function(e) {
    p
  })
  p + xlab("") + ylab("") + labs(color = "", fill = "") + 
    scale_alpha(guide = "none") + theme_light() + theme(legend.direction = "horizontal", 
                                                        legend.position = "top")
}

#' Plots unit by unit difference-in-differences. Dot size indicates the weights omega_i
#' used in the average that yields our treatment effect estimate. 
#' This estimate and endpoints of a 95% CI are plotted as horizontal lines.
#' Requires ggplot2
#' @param estimates as output by synthdid_estimate. Can be a single one or a list of them.
#' @param negligible.threshold Unit weight threshold below which units are plotted as small, transparent xs instead of circles. Defaults to .001.
#' @param negligible.alpha Determines transparency of those xs.
#' @param se.method the method used to calculate standard errors for the CI. See vcov.synthdid_estimate. 
#'        Defaults to 'jackknife' for speed. If 'none', don't plot a CI.
#' @param units a list of control units --- elements of rownames(Y) --- to plot differences for. Defaults to NULL, meaning all of them.
#' @export synthdid_units_plot
synthdid_units_plot = function(estimates, negligible.threshold = .001, negligible.alpha = .3, se.method='jackknife', 
                               se.manual=NULL, ci.manual=NULL, units=NULL) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .ignore <- tryCatch(attachNamespace("ggplot2"), error = function(e) e)
  } else {
    stop("Plotting requires the package `ggplot2`. Install it to use this function.")
  }
  if (class(estimates) == 'synthdid_estimate') { estimates = list(estimates) }
  if (is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
  plot.data = do.call(rbind, lapply(1:length(estimates), function(ee) {
    estimate = estimates[[ee]]
    setup = attr(estimate, 'setup')
    weights = attr(estimate, 'weights')
    Y = setup$Y - contract3(setup$X, weights$beta)
    N0 = setup$N0; N1 = nrow(Y) - N0
    T0 = setup$T0; T1 = ncol(Y) - T0
    
    lambda.pre = c(weights$lambda, rep(0, T1))
    lambda.post = c(rep(0, T0), rep(1 / T1, T1))
    omega.control = c(weights$omega, rep(0, N1))
    omega.treat = c(rep(0, N0), rep(1 / N1, N1))
    difs = as.vector(t(omega.treat) %*% Y %*% (lambda.post - lambda.pre)) - as.vector(Y[1:N0, ] %*% (lambda.post - lambda.pre))
    se = if (se.method == 'none') { NA } else { if (!is.null(se.manual)){
      se.manual[[ee]]
    } else{
      sqrt(vcov(estimate, method = se.method))
    }
      }
    include.units = if(is.null(units)) { 1:N0 } else { which(rownames(Y)[1:N0] %in% units) }
    data.frame(y = difs[include.units], unit = rownames(Y)[include.units], weight = omega.control[include.units],
               estimate = c(estimate), se = se, estimator = names(estimates)[[ee]], conf.low = ifelse(!is.null(ci.manual),ci.manual[ee*2-1],NA), conf.high = ifelse(!is.null(ci.manual),ci.manual[ee*2],NA))
  }))
  p = ggplot(plot.data) +
    geom_point(aes(x = unit, y = y, size = weight), data = plot.data[plot.data$weight > negligible.threshold, ]) +
    geom_point(aes(x = unit, y = y, size = weight), data = plot.data[plot.data$weight <= negligible.threshold, ], alpha = negligible.alpha, shape = 4, show.legend = FALSE) +
    geom_hline(aes(yintercept = estimate), size = .75)
  if (!all(is.na(plot.data$se))) {
    p = p + geom_hline(aes(yintercept = estimate - 1.96 * se), size = .5, alpha = .5) +
      geom_hline(aes(yintercept = estimate + 1.96 * se), size = .5, alpha = .5)
  } else if(!is.null(ci.manual)){
    p = p + geom_hline(aes(yintercept = conf.low), size = .5, alpha = .5) +
      geom_hline(aes(yintercept = conf.high), size = .5, alpha = .5)
  }
  p + facet_grid(. ~ estimator) + xlab('') + ylab('') + guides(shape = 'none') +
    theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


