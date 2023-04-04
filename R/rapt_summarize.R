#' rescale_pattern Function
#'
#' @description
#' Rescales an input 3D point patterns
#'
#' @details
#' Takes an input 3D point pattern (pp3 object) and scales it to have intensity of new.intensity by multiplying
#' the coordinates of each point as well as the window.
#'
#' @param pattern 3D point pattern (pp3 object) to be scaled
#' @param new_intensity final intensity of point pattern
#'
#' @return pp3_new will preserve all properties of the input \emph{pattern}, except the intensity
#' @export
rescale_pattern <- function(pattern, new_intensity) {
  intensity <- sum(spatstat.geom::intensity(pattern))
  if (is.pp3(pattern)) {
    factor <- (intensity / new_intensity)^(1 / 3)
    coo <- coords(pattern)
    coo_scaled <- coo * factor

    win_scaled <- box3(
      (domain(pattern)$xrange * factor),
      (domain(pattern)$yrange * factor),
      (domain(pattern)$zrange * factor)
    )
    out <- pp3(coo_scaled$x, coo_scaled$y, coo_scaled$z, win_scaled)
  } else if (is.ppp(pattern)) {
    factor <- (intensity / new_intensity)^(1 / 2)
    coo <- coords(pattern)
    coo_scaled <- coo * factor
    win_scaled <- owin(
      (domain(pattern)$xrange * factor),
      (domain(pattern)$yrange * factor)
    )
    out <- ppp(coo_scaled$x, coo_scaled$y, win_scaled)
  }
  marks(out) <- marks(pattern)
  return(out)
}


#' Average relabelings
#' @param relabelings output of  \code{\link[rTEM]{relabel_summarize}} function
#' @param envelope.value size of envelope to compute.  Should be decimal (e.g. 0.95 = 95\%)
#' @param funcs vector of summary functions to calculate
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description
#' Function take all relabelings and return averages and envelope values
#'
#' @details Take output of  \code{\link[rTEM]{relabel_summarize}} and create envelopes
#'
#' @return data frame with x value (distance), average y value, and y value envelopess
#' @export
average_relabelings <- function(relabelings, envelope.value = .95,
                                funcs = c("K", "G"),
                                K_cor = "trans", G_cor = "km", F_cor = "km",
                                GXGH_cor = "km", GXHG_cor = "km", G2_cor = "km") {
  # transform envelope value to high index (0.95 envelope will be 0.025 and 0.975)
  envelope.value <- envelope.value + (1 - envelope.value) / 2

  for (i in 2:10) {
    nm =paste("G", i, sep = "")
    if (any(funcs == nm)) {
      assign(paste("G", i, "_cor", sep=""), G2_cor)
    }
  }

  # get index of high and low envelope values and find values at each index
  hi.ind <- round((length(relabelings) + 1) * envelope.value, 0)
  lo.ind <- round((length(relabelings) + 1) * (1 - envelope.value), 0)
  if (lo.ind == 0) {
    lo.ind <- 1
  }
  # make the relabelings their own individual objects

  out = lapply(funcs, function(func) {
    cor = paste(func, "_cor", sep = "")
    cor = get(cor)

    rrl = paste("rrl_", func, sep = "")


    mmean <- sapply(relabelings, function(x) {
      x[[rrl]][[cor]]
    })

    ordered <- apply(mmean, 1, sort)

    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]

    # get r values
    r <- relabelings[[1]][[rrl]]$r
    # find the median at every distance
    med <- apply(mmean, 1, median)
    return(data.frame(r = r, mmean = med, lo = lo, hi = hi))
  })
  # K
  names(out) = sapply(1:length(funcs), function(func) {
    paste("rrl_", funcs[func], sep = "")
  })
  return(out)
}








#' Relabel input point pattern and calcluate summary functions
#' @param seed number to use in `set.seed` for reproducibility
#' @param funcs vector of summary functions to calculate
#' @param pattern point pattern of type ppp or pp3
#' @param dopant_formula mark of points from which to calculate summary functions
#' @param host_formula marks of points to not use in summary functions (except cross type functions)
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#' @param maxKr a numeric
#' @param nKr a numeric
#' @param maxGr a numeric
#' @param nGr a numeric
#' @param maxGXGHr a numeric
#' @param maxGXHGr a numeric
#' @param nGXr a numeric
#'
#' @description takes an input binary point pattern, randomly relabels it while maintainig the
#' same ratio of marks, and calculates summary functions.  Can then use  \code{\link[rTEM]{average_relabelings}}
#' to find averages and envelopes
#' @return summary functions listed in `funcs` variable
#' @export
relabel_summarize <- function(seed, pattern, funcs = c("K", "G", "F", "GXGH"),
                              dopant_formula, host_formula, k = 1,
                              maxKr = 10, nKr = 200, maxGr = 5, nGr = 1000,
                              maxGXGHr = 3, maxGXHGr = 8, nGXr = 1000, vside = 0.3,
                              K_cor = "trans", G_cor = "km", F_cor = "km",
                              GXGH_cor = "km", GXHG_cor = "km",
                              G2_cor = "km") {
  # set seed to make sure each random relabeling is different
  set.seed(seed)

  # calculate total host and dopant
  dopant_total <- sum(marks(pattern) == dopant_formula)
  host_total <- sum(marks(pattern) == host_formula)

  # relabel the input pattern, while maintaining the same proportion of host and dopant points
  relabeled <- rlabel(pattern,
                      labels = as.factor(c(
                        rep(dopant_formula, dopant_total),
                        rep(host_formula, host_total)
                      )), permute = TRUE
  )

  # select only the dopant type points
  dopant_relabeled <- subset(relabeled, marks == dopant_formula)

  # calculate summary functions
  if (is.pp3(pattern)) {
    if ("K" %in% funcs) {
      if (K_cor == "border") {
        rrl_K <- bK3est(dopant_relabeled, rmax = maxKr, nrval = nKr)
      } else {
        rrl_K <- K3est(dopant_relabeled, rmax = maxKr, nrval = nKr, correction = K_cor)
      }
    }
    if ("G" %in% funcs) {
      rrl_G <- G3est(dopant_relabeled, rmax = maxGr, nrval = nGr, correction = G_cor)
    }

    for (i in 2:10) {
      if (paste("G", i, sep = "") %in% funcs) {
        G_i = G3est_nn(dopant_relabeled,  rmax = maxGr, nrval = nGr,
                       k = i, correction = G2_cor)
        assign(paste("rrl_G", i, sep=""), G_i)
      }
    }

    if ("F" %in% funcs) {
      rrl_F <- F3est(dopant_relabeled, rmax = maxGr, nrval = nGr, correction = F_cor, vside = vside)
    }
    if ("GXGH" %in% funcs) {
      rrl_GXGH <- G3cross(relabeled,
                          i = dopant_formula, j = host_formula,
                          rmax = maxGXGHr, nrval = nGXr, correction = GXGH_cor
      )
    }
    if ("GXHG" %in% funcs) {
      rrl_GXHG <- G3cross(relabeled,
                          i = host_formula, j = dopant_formula,
                          rmax = maxGXHGr, nrval = nGXr, correction = GXHG_cor
      )
    }
  } else if (is.ppp(pattern)) {
    if ("K" %in% funcs) {
      rrl_K <- Kest(dopant_relabeled, r = seq(0, maxKr, length.out = nKr), correction = K_cor)
    }
    if ("G" %in% funcs) {
      rrl_G <- Gest_nn(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = G_cor, k = 1)

    }
    for (i in 2:10) {
      if (paste("G", i, sep = "") %in% funcs) {
        G_i = Gest_nn(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = G2_cor, k = i)
        assign(paste("rrl_G", i, sep=""), G_i)
      }
    }
    if ("F" %in% funcs) {
      rrl_F <- Fest(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = F_cor)
    }
    if ("GXGH" %in% funcs) {
      rrl_GXGH <- Gcross(relabeled,
                         i = dopant_formula, j = host_formula,
                         r = seq(0, maxGXGHr, length.out = nGXr), correction = GXGH_cor
      )
    }
    if ("GXHG" %in% funcs) {
      rrl_GXHG <- Gcross(relabeled,
                         i = host_formula, j = dopant_formula,
                         r = seq(0, maxGXHGr, length.out = nGXr), correction = GXHG_cor
      )
    }
  }


  all_funcs <- c("K", "G", "F", "GXGH", "GXHG")

  G_nn = rep(NA, 9)
  for (i in 2:10) {
    G_nn[i-1] = paste("G", i, sep = "")
  }
  all_funcs = c(all_funcs, G_nn)

  #all_funcs %in% funcs
  relabs <- c("rrl_K", "rrl_G", "rrl_F", "rrl_GXGH", "rrl_GXHG")
  relabs = c(relabs, G_nn)
  out <- lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) <- relabs[all_funcs %in% funcs]
  return(out)
}


#' function to calculate all summary functions on point pattern
#' @param funcs vector of summary functions to calculate
#' @param pattern point pattern of type ppp or pp3
#' @param dopant_formula mark of points from which to calculate summary functions
#' @param host_formula marks of points to not use in summary functions (except cross type functions)
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#' @param ... maxKr, nKr, maxGr, nGr, maxGXGHr, maxGXHGr, nGXr
#'
#' @description Cacluate that summary functions included in `funcs` on `pattern`, a point pattern
#' of class ppp or pp3.  Now can do advanced nearest neighbor if `funcs` contains `G2`, `G3`, etc.
#' @export
calc_summary_funcs <- function(pattern, funcs = c("K", "G", "F", "GXGH"),
                               dopant_formula, host_formula,
                               maxKr = 10, nKr = 200, maxGr = 5, nGr = 1000,
                               maxGXGHr = 3, maxGXHGr = 8, nGXr = 1000, vside = 0.3,
                               K_cor = "trans", G_cor = "km", F_cor = "km",
                               GXGH_cor = "km", GXHG_cor = "km", G2_cor = "km") {
  # select only the dopant type points
  pattern_dopant <- subset(pattern, marks == dopant_formula)

  if (is.pp3(pattern)) {
    # calculate summary functions
    if ("K" %in% funcs) {
      if (K_cor == "border") {
        K <- bK3est(pattern_dopant, rmax = maxKr, nrval = nKr)
      } else {
        K <- K3est(pattern_dopant, rmax = maxKr, nrval = nKr, correction = K_cor)
      }
    }
    if ("G" %in% funcs) {
      G <- G3est(pattern_dopant, rmax = maxGr, nrval = nGr, correction = G_cor)
    }
    for (i in 2:10) {
      if (paste("G", i, sep = "") %in% funcs) {
        G_i = G3est_nn(dopant_relabeled,  rmax = maxGr, nrval = nGr,
                       k = i, correction = G2_cor)
        assign(paste("G", i, sep=""), G_i)
      }
    }
    if ("F" %in% funcs) {
      F <- F3est(pattern_dopant, rmax = maxGr, nrval = nGr, correction = F_cor, vside = vside)
    }
    if ("GXGH" %in% funcs) {
      GXGH <- G3cross(pattern,
                      i = dopant_formula, j = host_formula,
                      rmax = maxGXGHr, nrval = nGXr, correction = GXGH_cor
      )
    }
    if ("GXHG" %in% funcs) {
      GXHG <- G3cross(pattern,
                      i = host_formula, j = dopant_formula,
                      rmax = maxGXHGr, nrval = nGXr, correction = GXHG_cor
      )
    }
  } else if (is.ppp(pattern)) {
    if ("K" %in% funcs) {
      K <- Kest(pattern_dopant, r = seq(0, maxKr, length.out = nKr), correction = K_cor)
    }
    if ("G" %in% funcs) {
      G <- Gest_nn(pattern_dopant, r = seq(0, maxGr, length.out = nGr), correction = G_cor, k = 1)
    }
    for (i in 2:10) {
      if (paste("G", i, sep = "") %in% funcs) {
        G_i = Gest_nn(pattern_dopant, r = seq(0, maxGr, length.out = nGr),
                      correction = G2_cor, k = i)
        assign(paste("G", i, sep=""), G_i)
      }
    }
    if ("F" %in% funcs) {
      F <- Fest(pattern_dopant, r = seq(0, maxGr, length.out = nGr), correction = F_cor)
    }
    if ("GXGH" %in% funcs) {
      GXGH <- Gcross(pattern,
                     i = dopant_formula, j = host_formula,
                     r = seq(0, maxGXGHr, length.out = nGXr), correction = GXGH_cor
      )
    }
    if ("GXHG" %in% funcs) {
      GXHG <- Gcross(pattern,
                     i = host_formula, j = dopant_formula,
                     r = seq(0, maxGXHGr, length.out = nGXr), correction = GXHG_cor
      )
    }
  }


  all_funcs <- c("K", "G", "F", "GXGH", "GXHG")
  G_nn = rep(NA, 9)
  for (i in 2:10) {
    G_nn[i-1] = paste("G", i, sep = "")
  }
  all_funcs = c(all_funcs, G_nn)
  relabs <- c("K", "G", "F", "GXGH", "GXHG")
  relabs = c(relabs, G_nn)
  out <- lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) <- relabs[all_funcs %in% funcs]
  return(out)
}



#' Plot summary functions (developement version)
#' @param func Summary function to plot.
#' @param observed_values a list. observed summary function value. Must have one of the following structures
#' \itemize{
#'  \item{Option 1: List of summary functions}{observed_values[[func]][[cor]]}
#'  \item{Option 2: List of lists of summary functions}{observed_values[[pattern_num]][[func]][[cor]]}
#'  \item{Option 2: Dataframe}{observed_values[,c("r", y_value)]}
#'  }
#' @param envelopes a list envelope values found by
#' function such as   \code{\link{calc_summary_funcs}}.  Should be list of length equal
#' to the number of envelopes.  The structure should resemble:
#' `envelopes[[envelope_num]]$rrl_K[, c("r", "mmean")]`
#' @param pattern.colors colors and names for envelope lines.
#'  MUST follow some formatting rules.  Envelope names must come first.
#'  If each observed value corresponds to a different envelope, then
#'  the number of `observed_values` must match `length(envelopes)` and only
#'  the first `1:length(envelopes)` of `pattern_colors` will be used.
#'  and must set `base_value = "each"`.  If not, each envelope and each observed
#'  value needs its own color.  then the `mmean` value from
#'  `envelopes[[1]]` will be subtracted from each envelope and observed value.
#' @param fill.colors colors and names for envelope fill.  Match names to `pattern.colors`.
#' Recommended to leave as NA and it will automattically be set to match  `pattern.colors`
#' @param sqrt either `"K", "all", or "none"`. If `"K"`,
#' then only \code{expression(tilde(K)[g](r))} will be found using the differences of the
#' square roots, rather than differences.  If `"all"`, then all functions will be found
#' using such.  If `"none"` (or actually anything other than the first two options) then
#' no square roots are taken. Ignored if `raw = "TRUE"`
#' @param raw.  A logical.  If TRUE, then no envelope mmeans are subtracted from each value
#' @param base_value a character, either `"first" or "each"`.  Does each observed value
#' correspond to a different envelope?  If yes, set to `"each"`. Otherwise, set to
#' `"first"` and the envelopes[[1]] will be used for each
#' @param unit a character.  This will appear as the units in the x axis label
#' @param K_cor edge correction(s) to be used when `func = "K"`
#' @param G_cor edge correction(s) to be used when `func = "G"`
#' @param F_cor edge correction(s) to be used when `func = "F"`
#' @param GXGH_cor edge correction(s) to be used when `func = "GXGH"`
#' @param GXHG_cor edge correction(s) to be used when `func = "GXHG"`
#' @param alpha numeric between 0 and 1.  Transparency of envelopes
#' @param legend.key.size numeric.  size of legend key
#' @param legend.text.size numeric. size of legend text
#' @param legend.position.  vector of 2 numerics.  coordinates of legend
#' @param axis.title.size numeric.  Size of axis title
#' @param title a character.  Text for title
#' @param title.size numeric.  Title size
#' @param axis.text.x.size numeric.  size of text on x axis
#' @param axis.text.y.size numeric.  size of text on y axis
#' @param linewidth numeric.  Width of lines in plot
#' @param env_linewidth numeric.  Width of lines that make up envelope edges
#' @param linetype a character.  Type of lines that make up lines
#' @param env_linetype a character.  Type of lines that make up  envelope lines
#' @description Plot the observed value with envelopes for expected values for summary function
#' @details The best way to learn about this function is to read the parameter definitions.
#' @export
plot_summary <- function(func = "K",
                         observed_values, envelopes, ...,
                         pattern.colors = c("Envelope 1" = "pink", "Envelope 2" = "gray", "Observed" = "blue"),
                         fill.colors =  NA,
                         sqrt = "K", raw = "FALSE",
                         base_value = "first", unit = "nm",
                         K_cor = "trans", G_cor = "km", F_cor = "km",
                         GXGH_cor = "km", GXHG_cor = "km", G2_cor = "km",
                         alpha = 0.5,
                         legend.key.size = 20, legend.text.size = 20,
                         legend.position =  c(0.75, 0.80),
                         axis.title.size = 40,
                         title = "", title.size = 40, axis.text.x.size = 40,
                         axis.text.y.size = 40, linewidth = 0.5,  env_linewidth = 0.5,
                         linetype = "solid",
                         env_linetype = "dashed") {
  if (all(is.na(fill.colors))) {
    fill.colors = pattern.colors
  }


  for (i in 2:10) {
    nm =paste("G", i, sep = "")
    if (func == nm) {
      assign(paste("G", i, "_cor", sep=""), G2_cor)
    }
  }

  cor = paste(func, "_cor", sep = "")

  cor = get(cor)

  rrl = paste("rrl_", func, sep = "")
  ylabs = list("K" = expression(tilde(K)[g](r)), "G" = expression(tilde(G)[g](r)),
               "F" = expression(tilde(F)[g](r)),
               "GXGH" = expression(tilde(G)[gh](r)), "GXHG" = expression(tilde(G)[hg](r)),
               "G2" = expression(tilde(G)[2,g](r)), "G3" = expression(tilde(G)[3,g](r)),
               "G4" = expression(tilde(G)[4,g](r)), "G5" = expression(tilde(G)[5,g](r)),
               "G6" = expression(tilde(G)[6,g](r)), "G7" = expression(tilde(G)[7,g](r)),
               "G8" = expression(tilde(G)[8,g](r)), "G9" = expression(tilde(G)[9,g](r)),
               "G10" = expression(tilde(G)[10,g](r)))


  ## if returning the square root of the difference, rather than the difference
  if ((sqrt == "all") || (sqrt == "K" && func == "K")) {
    take_root = function(vector) {
      return(sqrt(vector))
    }
  }

  # if not, just return the original vector
  else {
    take_root = function(vector) {
      return(vector)
    }
  }
  # if there is an envelope for each observed value, then plot envelopes separately
  if ((length(envelopes) == length(observed_values)) && base_value != "first") {
    print("start each")
    long = data.frame("r" = c(),
                      "mmean" = c(),
                      "lo" = c(),
                      "hi" = c(),
                      "type" = c())
    i <- 1
    while (i <= length(envelopes)) {
      temp <- envelopes[[i]][[rrl]]
      observed <- observed_values[[i]][[func]]
      temp$type <- names(pattern.colors)[i]
      baseline = take_root(temp$mmean)

      if (raw) {
        baseline = 0
      }

      temp$lo = take_root(temp$lo) - baseline
      temp$hi = take_root(temp$hi) - baseline
      temp$mmean = take_root(observed[[cor]]) - baseline

      long <- rbind(long, temp)
      i <- i + 1
    }

    gplot <- long %>% ggplot2::ggplot(aes(
      x = r, ymin = lo,
      ymax = hi, color = type, fill = type
    )) +
      geom_ribbon(alpha = alpha, linewidth = env_linewidth, linetype = env_linetype) +
      geom_hline(yintercept = 0) +
      geom_line(aes(x= r, y = mmean, color = type), linewidth = linewidth, linetype = linetype)


    gplot <- gplot + theme(
      plot.title = element_text(hjust = 0.5, size = title.size),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = axis.text.x.size),
      axis.text.y = element_text(size = axis.text.y.size),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = axis.title.size),
      legend.key.size = unit(legend.key.size, "pt"),
      legend.text = element_text(size = legend.text.size),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position =legend.position,
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),# ...
    ) +
      # guides(color = "none") +
      scale_color_manual(values = pattern.colors) +scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(ylabs[[func]]) +
      ggtitle(title)


  }

  #######

  else {
    print("start first")
    baseline = envelopes[[1]][[rrl]]$mmean
    if (raw) {
      baseline = 0
    }
    long = data.frame("r" = c(),
                      "mmean" = c(),
                      "lo" = c(),
                      "hi" = c(),
                      "type" = c())

    i <- 1
    while (i <= length(envelopes)) {
      temp <- envelopes[[i]][[rrl]]
      temp$type <- names(pattern.colors)[i]
      temp$lo = take_root(temp$lo) - take_root(baseline)
      temp$hi = take_root(temp$hi) - take_root(baseline)
      temp$mmean = take_root(temp$mmean) - take_root(baseline)
      long <- rbind(long, temp)
      i <- i + 1

    }

    #class(observed_funcs[[image_num]][[func]])
    if (is.data.frame(observed_values)) {
      observed <- data.frame(
        r = observed_values[["r"]],

        mmean = take_root(observed_values[[cor]]) - take_root(baseline),
        lo = take_root(observed_values[[cor]]) - take_root(baseline),
        hi = take_root(observed_values[[cor]]) - take_root(baseline),
        type = "Observed"
      )
    }


    else if (is.list(observed_values)) {
      if (is.data.frame(observed_values[[1]])) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        obs = observed_values[[func]]
        obs_01 <- data.frame(
          r = obs$r,

          mmean = take_root(obs[[cor]]) - take_root(baseline),
          lo = take_root(obs[[cor]]) - take_root(baseline),
          hi = take_root(obs[[cor]]) - take_root(baseline),
          type = names(pattern.colors)[length(envelopes) + 1])
        observed = rbind(observed, obs_01)
      }

      else if (is.list(observed_values[[1]])) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        for (i in 1:length(observed_values)) { #
          obs = (observed_values[[i]][[func]])
          obs_01 <- data.frame(
            r = obs$r,

            mmean = take_root(obs[[cor]]) - take_root(baseline),
            lo = take_root(obs[[cor]]) - take_root(baseline),
            hi = take_root(obs[[cor]]) - take_root(baseline),
            type = names(pattern.colors)[length(envelopes) + i])
          observed = rbind(observed, obs_01)
        }
      }


    }

    else {
      stop("observed_values must be either dataframe or list of dataframes with a single summary function,
             list of dataframes, each dataframe containing the same function,
             or a list of list of dataframes, the outer list being the pattern, the inner list being the different summary functions")
    }


    long <- rbind(long, observed)

    gplot <- long %>% ggplot2::ggplot(aes(
      x = r, ymin = (lo) ,
      ymax = (hi) , color = type, fill = type
    )) +
      geom_ribbon(alpha = alpha, linewidth = linewidth, linetype = linetype) +
      geom_hline(yintercept = 0)# +
    # geom_line(aes(x= r, y = sqrt(mmean) , color = type))


    print("start plot")
    gplot <- gplot + theme(
      plot.title = element_text(hjust = 0.5, size = title.size),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = axis.text.x.size),
      axis.text.y = element_text(size = axis.text.y.size),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = axis.title.size),
      legend.key.size = unit(legend.key.size, "pt"),
      legend.text = element_text(size = legend.text.size),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position =legend.position,
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),# ...
    ) +
      # guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(ylabs[[func]]) +
      ggtitle(title)
  }
}

