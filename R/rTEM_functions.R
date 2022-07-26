#' Rescaler Function
#'
#'@description
#' Rescales an input 3D point pattern
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
rescaler = function(pattern, new_intensity) {
  intensity = spatstat.geom::intensity(pattern)
  factor = (intensity / new_intensity)^(1/3)
  coo = coords(pattern)
  coo_scaled = coo*factor


  win_scaled = box3((domain(pattern)$xrange*factor),
                    (domain(pattern)$yrange*factor),
                    (domain(pattern)$zrange*factor))
  pp3_new = pp3(coo_scaled$x, coo_scaled$y, coo_scaled$z, win_scaled)
}

#' Calculate T value
#'
#'
#' @description
#' Takes input of relabeling object and calculates the T value for either K or G function
#'
#' @details This uses the method from \emph{Baddeley et al.} to calculate the T value for an envelope for either
#' K  (\emph{K3est}) G \emph{G3est} function
#' @references
#' Baddeley, A., Diggle, P. J., Hardegen, A., Lawrence, T_, Milne, R. K., & Nair, G. (2014).
#' On tests of spatial pattern based on simulation envelopes. Ecological Monographs, 84(3), 477–489. https://doi.org/10.1890/13-2042.1
#' @export
Tcalc = function(relabelings, func = "K", rmin = 0, rmax, K_cor = "border", G_cor = "km") {
  if (func == "K") {
    if (K_cor == "trans") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_K$trans
      })
    }
    else if (K_cor == "iso") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_K$iso
      })
    }
    else if (K_cor == "border") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_K$bord
      })
    }
    else {
      print("Incorrect K edge correction")
    }

    r = relabelings[[1]]$rrl_K$r
    med = apply(mmean, 1, median)
    ind = which(r <=  rmax & r >= rmin)
    interval = rmax / (length(r) -1)
    T = apply(mmean, 2, function(x) {
      T_1= sqrt(x) - sqrt(med)
      T_1 = T_1[ind]
      sum(T_1^2) * interval
    })
    return(T)
  }
  if (func == "G") {
    if (G_cor == "km") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_G$km
      })
    }

    else if (G_cor == "rs") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_G$rs
      })
    }
    else {
      print("Incorrect G edge correction")
    }
    r = relabelings[[1]]$rrl_G$r
    med = apply(mmean, 1, median)
    ind = which(r <=  rmax & r >= rmin)
    interval = rmax / (length(r) -1)
    T = apply(mmean, 2, function(x) {
      T_1=x - med
      T_1 = T_1[ind]
      sum(T_1^2) * interval
    })
  }
  return(T)
}

#' Functionto take all relabelings and return averages and envelope values
#' @param relabelings output of `relabeler` function
#' @param envelope.value
#' @export
rrl_averager = function(relabelings, envelope.value = .95,
                        funcs = c("K", "G"),
                        K_cor = "trans", G_cor = "km", F_cor = "km",
                        GXGH_cor = "km", GXHG_cor = "km") {

  # transform envelope value to high index (0.95 envelope will be 0.025 and 0.975)
  envelope.value =  envelope.value + (1- envelope.value)/2

  # make the relabelings their own individual objects

  # K
  # extract K(r) values
  if ("K" %in% funcs) {
    if (K_cor == "trans") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_K$trans
      })
    }
    else if (K_cor == "iso") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_K$iso
      })
    }
    else if (K_cor == "border") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_K$bord
      })
    }
    else {
      print("Incorrect K edge correction")
    }

    # order K(r) values by value
    ordered = apply(mmean, 1, sort)

    # get index of high and low envelope values and find values at each index
    hi.ind = round(length(relabelings) * envelope.value, 0)
    lo.ind = round(length(relabelings) * (1-envelope.value), 0)
    if (lo.ind == 0) {
      lo.ind = 1
    }
    lo = ordered[lo.ind,]
    hi = ordered[hi.ind,]

    # get r values
    r = relabelings[[1]]$rrl_K$r

    # find the median at every distance
    med = apply(mmean, 1, median)
    rrl_K = data.frame(r = r, mmean =med, lo = lo, hi = hi)

  }

  # Repeat for G, GX, and F
  # G
  if ("G" %in% funcs) {
    if (G_cor == "km") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_G$km
      })
    }

    else if (G_cor == "rs") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_G$rs
      })
    }
    else {
      print("Incorrect G edge correction")
    }

    r = relabelings[[1]]$rrl_G$r
    ordered = apply(mmean, 1, sort)
    lo = ordered[lo.ind,]
    hi = ordered[hi.ind,]
    med = apply(mmean, 1, median)
    rrl_G = data.frame(r = r, mmean =med, lo = lo, hi = hi)
  }

  # F
  if ("F" %in% funcs) {

    if (F_cor == "km") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_F$km
      })
    }

    else if (F_cor == "rs") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_F$rs
      })
    }
    else if (F_cor == "cs") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_F$cs
      })
    }

    else {
      print("Incorrect F edge correction")
    }
    r = relabelings[[1]]$rrl_F$r
    ordered = apply(mmean, 1, sort)
    lo = ordered[lo.ind,]
    hi = ordered[hi.ind,]
    med = apply(mmean, 1, median)
    rrl_F = data.frame(r = r, mmean =med, lo = lo, hi = hi)

  }

  # GXGH
  if ("GXGH" %in% funcs) {
    if (GXGH_cor == "km") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_GXGH$km
      })
    }

    else if (GXGH_cor == "rs") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_GXGH$rs
      })
    }
    else if (GXGH_cor == "han") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_GXGH$han
      })
    }
    else if (GXGH_cor == "none") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_GXGH$raw
      })
    }

    else {
      print("Incorrect GXGH edge correction")
    }

    r = relabelings[[1]]$rrl_GXGH$r
    ordered = apply(mmean, 1, sort)
    lo = ordered[lo.ind,]
    hi = ordered[hi.ind,]
    med = apply(mmean, 1, median)
    rrl_GXGH = data.frame(r = r, mmean =med, lo = lo, hi = hi)
  }

  # GXGH
  if ("GXHG" %in% funcs) {
    if (GXHG_cor == "km") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_GXHG$km
      })
    }

    else if (GXHG_cor == "rs") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_GXHG$rs
      })
    }
    else if (GXHG_cor == "han") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_GXHG$han
      })
    }
    else if (GXHG_cor == "none") {
      mmean =sapply(relabelings, function(x) {
        x$rrl_GXHG$raw
      })
    }

    else {
      print("Incorrect GXHG edge correction")
    }

    r = relabelings[[1]]$rrl_GXHG$r
    ordered = apply(mmean, 1, sort)
    lo = ordered[lo.ind,]
    hi = ordered[hi.ind,]
    med = apply(mmean, 1, median)
    rrl_GXHG = data.frame(r = r, mmean =med, lo = lo, hi = hi)

  }
  all_funcs = c("K", "G", "F", "GXGH", "GXHG")
  all_funcs %in% funcs
  relabs = c("rrl_K", "rrl_G", "rrl_F", "rrl_GXGH", "rrl_GXHG")
  out =lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) =  relabs[all_funcs %in% funcs]
  return(out)
}




#' Simulate 2D Multimers
#' @description Simulate
multimersim = function(exp_ppp,  thickness, group_size = 2,
                       num_neighbors = 6, neighbor_number = 1, weights = c(1, 1, 1/4), probs = c(0.6, 0.2, 0.2, 0),
                       rcp_pattern, intensity_rcp = 1, maxGr, maxKr, nGr, nKr) {
  ## make probability vector same length as num_neighbors
  if (num_neighbors > length(probs)) {
    probs = c(probs, rep(0, num_neighbors - length(probs)))
  }

  # get desired intensity
  npoints = exp_ppp$n
  print(npoints)
  box_2d = exp_ppp$window

  box_area = area(box_2d)
  #vol = box_area * thickness
  # intensity_exp = npoints / vol
  # rescale RCP pattern to match physical system (1 point per nm)
  rcp_scaled = rescaler(rcp_pattern, intensity_rcp)

  # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
  rcp_box = subset(rcp_scaled, x > box_2d$xrange[1] & x < box_2d$xrange[2] &
                     y > box_2d$yrange[1] & y < box_2d$yrange[2] &
                     z >0 & z < thickness)
  ## label n/2 points as dopant and the rest as host
  n_irp = npoints
  n_total = length(rcp_box$data$x)
  n_host = n_total - n_irp

  # determine the number of groups of molecules to be present (if multimers, then n_irp/2)
  n_groups = round(n_irp/ group_size,0)
  rcp_labeled = rlabel(rcp_box, labels = c(rep("A",n_groups) , rep("C", length(rcp_box$data$x) - n_groups)), permute = TRUE)
  # extract dopant points
  rcp_dope = subset(rcp_labeled, marks == "A")
  # extract host points
  rcp_host = subset(rcp_labeled, marks == "C")
  # find 6 nearest host to each dopant
  nn_ind= lapply(1:num_neighbors, function(x) {
    nncross(rcp_dope, rcp_host, k = x)
  })

  # get nn distance x, y, and z values for each neighbor
  dist_coords = lapply(1:num_neighbors, function(x) {
    coords(rcp_dope) - coords(rcp_host[nn_ind[[x]][,2]])
  })
  # this just rearranges dist_coords to be grouped by point rather than neighbor order
  holder = c()
  neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
    for (i in 1:length(dist_coords)) {
      holder = rbind(holder, dist_coords[[i]][x,])
    }
    return(holder)
  })

  # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
  # rank will be from smallest to largest
  ranks = lapply(1:length(neighbors), function(outer) {
    distances = sapply(1:nrow(neighbors[[outer]]), function(inner) {
      sqrt(weights[1]*neighbors[[outer]]$x[inner]^2 + weights[2]*neighbors[[outer]]$y[inner]^2 +weights[3]*neighbors[[outer]]$z[inner]^2)
    })
    order(distances)
    rank(distances)
  })

  # semi randomly select one of the neighbors, prefering onces with lower x-y distance
  which_neighbor = t(sapply(1:length(ranks), function(i) {
    sample(1:num_neighbors, size =group_size-1, prob = probs[ranks[[i]]])
  }))
  if (group_size == 2) {
    which_neighbor = which_neighbor[1,]
  }

  which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)

  # get the coordinates of the semi randomly selected neighbor
  ind= t(sapply(1:nrow(which_neighbor), function(i) {
    sapply(2:(group_size), function(j) {
      nn_ind[[which_neighbor[i,j]]][i,2]
    })
  }))

  #duplicate = apply(ind, 2, duplicated)
  #duplicate = duplicated(c(ind))
  duplicate = duplicated(ind, MARGIN = 0)
  duplicate_ind = which(duplicate, arr.ind = TRUE)[,1]
  duplicate_ind = unique(duplicate_ind)
  not_duplicate = ind[!duplicate]
  ## points that are not duplicates
  chosen_points = rcp_host$data[not_duplicate]


  ## host pattern with previously used points removed
  host_unique_pp3 = rcp_host[-ind]
  #unique_points = rcp_host$data[-ind[duplicate],]
  print(paste("first duplicate ", sum(duplicate)))


  ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two dopants )
  # take the 2nd nearest neighbor
  while (sum(duplicate >0) ) {

    # find 6 nearest host to each dopant
    next_nn= lapply(1:num_neighbors, function(x) {
      nncross(rcp_dope, host_unique_pp3, k = x)
    })

    ####
    dist_coords = lapply(1:num_neighbors, function(x) {
      coords(rcp_dope) - coords(host_unique_pp3[next_nn[[x]][,2]])
    })
    ####

    # this just rearranges dist_coords to be grouped by point rather than neighbor order
    holder = c()
    neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
      for (i in 1:length(dist_coords)) {
        holder = rbind(holder, dist_coords[[i]][x,])
      }
      return(holder)
    })

    # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
    # rank will be from smallest to largest
    ranks = lapply(1:length(neighbors), function(outer) {
      distances = sapply(1:nrow(neighbors[[outer]]), function(inner) {
        sqrt(weights[1]*neighbors[[outer]]$x[inner]^2 + weights[2]*neighbors[[outer]]$y[inner]^2 +weights[3]*neighbors[[outer]]$z[inner]^2)
      })
      order(distances)
      rank(distances)
    })


    # semi randomly select one of the neighbors, prefering onces with lower x-y distance
    which_neighbor = t(sapply(1:length(ranks), function(i) {
      sample(1:num_neighbors, size =group_size-1, prob = probs[ranks[[i]]])
    }))
    if (group_size == 2) {
      which_neighbor = which_neighbor[1,]
    }
    which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)


    # get the coordinates of the semi randomly selected neighbor
    all_ind= t(sapply(1:nrow(which_neighbor), function(i) {
      sapply(2:(group_size), function(j) {
        next_nn[[which_neighbor[i,j]]][i,2]
      })
    }))

    #ind[duplicate] = all_ind[duplicate]
    #all_ind = c(all_ind)
    # duplicate = c(duplicate)
    # extract just the ones needed to replace duplicates
    next_ind =all_ind[duplicate]

    ## use two duplicates - one that has the index of duplicates in entire all_ind, one that has actual new duplicates
    #new = all_ind[duplicate]
    mat_1 = matrix(FALSE, nrow = nrow(all_ind), ncol = ncol(all_ind))
    mat_1[duplicate] = all_ind[duplicate]
    # zeros = test ==0
    duplicate_mat = duplicated(mat_1, MARGIN = 0)
    duplicate_mat[mat_1==0] = FALSE


    ## get duplicates and their index
    duplicate = duplicated(next_ind, MARGIN = 0)
    duplicate_ind = which(duplicate, arr.ind = TRUE)
    duplicate_ind = unique(duplicate_ind)
    not_duplicate = unique(next_ind)
    ## points that are not duplicates
    new_points = host_unique_pp3$data[not_duplicate]
    chosen_points = rbind(chosen_points, new_points)

    ## host pattern with previously used points removed
    host_unique_pp3 = host_unique_pp3[-next_ind]

    ## change duplicate back to full matrix form
    duplicate = duplicate_mat

    print(paste("loop duplicate ", sum(duplicate)))
  }
  dim(chosen_points)

  # get coordinates of all original points and their selected neighbors
  coords_multimer = rbind(coords(rcp_dope), data.frame(x = chosen_points$x, y = chosen_points$y, z =chosen_points$z))
  dim(coords_multimer)

  # make ppp and caluclate summary functions
  rcp_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
  rcp_G = Gest_nn(rcp_box, correction = "km", k = neighbor_number, r = seq(0, maxGr, length.out = nGr))
  rcp_K =Kest(rcp_box, r = seq(0, maxKr, length.out= nKr))
  summary = list(rrl_G = rcp_G, rrl_K = rcp_K)
  return(summary)
}

#' Relabel multimers
#'
#' multimer
#' @export
multimers_relab = function (exp_ppp, num_relabs, thickness, rcp_list, maxGr, maxKr,
                            nGr, nKr,
                            group_size = 2, num_neighbors = 6, neighbor_number = 1, weights = c(1, 1, 1/4),
                            probs = c(0.6, 0.2, 0.2, 0),
                            ncores = detectCores()) {

  ## make vectors of maximum radii to take T tests to for K and G
  G_rad_vec = seq(0, maxGr, length.out = (nGr/50) +1)
  G_rad_vec = G_rad_vec[2:length(G_rad_vec)]
  K_rad_vec = seq(0, maxKr, length.out = (nKr/50)+1)
  K_rad_vec = K_rad_vec[2:length(K_rad_vec)]

  ## number of relabelings per rcp pattern
  nper = num_relabs/floor(length(rcp_list))
  # index for which rcp pattern to use
  ind = rep(1:length(rcp_list), nper)

  cl = makeForkCluster(ncores, outfile = "outfile_test")
  ##  calculate envelopes for length(ind) rcp multimer patterns
  multimers = parLapply(cl, 1:length(ind), function(x) {
    print(paste("start rep  ", x))
    val = multimersim(exp_ppp, thickness = thickness, group_size = group_size,
                      num_neighbors =num_neighbors, neighbor_number = neighbor_number, weights = weights, probs = probs,
                      rcp_list[[ind[x]]],
                      intensity_rcp = 1, maxGr = maxGr, maxKr = maxKr, nGr = nGr, nKr = nKr)
    print(paste("end rep  ", x))
    val
  })

  clusterExport(cl, "multimers", envir = environment())
  # calculate T values for K and G summaries for a variety of max radii
  multimers_T_G = parLapply(cl, G_rad_vec, function(x) {
    Tcalc(multimers, x, func = "G", rmax =maxGr)
  })
  multimers_T_K = parLapply(cl, K_rad_vec, function(x) {
    Tcalc(multimers, x, func = "K",rmax = maxKr)
  })

  # get the 95% CI envelope
  multimers = rrl_averager(multimers, envelope_value = 0.95,  K_cor = "border", G_cor = "km") # use same name to get rid of the older (and enormous) object
  stopCluster(cl)
  return(list(relabs = multimers, T_G = multimers_T_G, T_K = multimers_T_K))
}




## functions for RCP
#' @export
func_rcp_simple = function(exp_ppp, thickness, rcp_pattern, neighbor_number = 1,
                           intensity_rcp = 1, maxGr, maxKr, nGr, nKr) {
  # get desired intensity
  npoints = exp_ppp$n
  print(npoints)
  box_2d = exp_ppp$window
  box_area = area(box_2d)
  #vol = box_area * thickness
  # intensity_exp = npoints / vol
  # rescale RCP pattern to match physical system (1 point per nm)
  rcp_scaled = rescaler(rcp_pattern, intensity_rcp)

  # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
  rcp_box = subset(rcp_scaled, x > box_2d$xrange[1] & x < box_2d$xrange[2] &
                     y > box_2d$yrange[1] & y < box_2d$yrange[2] &
                     z >0 & z < thickness)
  ## label n/2 points as dopant and the rest as host
  n_irp = npoints
  n_total = length(rcp_box$data$x)
  n_host = n_total - n_irp
  rcp_labeled = rlabel(rcp_box, labels = c(rep("A", n_irp) , rep("C", n_host)), permute = TRUE)
  # extract dopant points
  rcp_dope = subset(rcp_labeled, marks == "A")
  # extract host points
  #rcp_host = subset(rcp_labeled, marks == "C")

  coords_dope = coords(rcp_dope)
  #win_box = box3(box_2d$xrange, box_2d$yrange, c(0, thickness))
  rcp_box = ppp(x = coords_dope$x, y = coords_dope$y, window = box_2d)

  rcp_G = Gest_nn(rcp_box, correction = "km", k =neighbor_number, r = seq(0, maxGr, length.out = nGr))
  rcp_K =Kest(rcp_box, r = seq(0, maxKr, length.out= nKr))
  summary = list(rrl_G = rcp_G, rrl_K = rcp_K)
  return(summary)
}

#' Create  RCP Pattern with same intensity and number of points as experimental pattern
#' @export
func_rcp_simple = function(exp_ppp, thickness, rcp_pattern, neighbor_number = 1,
                           intensity_rcp = 1, maxGr, maxKr, nGr, nKr) {
  # get desired intensity
  npoints = exp_ppp$n
  print(npoints)
  box_2d = exp_ppp$window
  box_area = area(box_2d)

  #vol = box_area * thickness
  # intensity_exp = npoints / vol
  # rescale RCP pattern to match physical system (1 point per nm)
  rcp_scaled = rescaler(rcp_pattern, intensity_rcp)

  # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
  rcp_box = subset(rcp_scaled, x > box_2d$xrange[1] & x < box_2d$xrange[2] &
                     y > box_2d$yrange[1] & y < box_2d$yrange[2] &
                     z >0 & z < thickness)
  ## label n/2 points as dopant and the rest as host
  n_irp = npoints
  n_total = length(rcp_box$data$x)
  n_host = n_total - n_irp
  rcp_labeled = rlabel(rcp_box, labels = c(rep("A", n_irp) , rep("C", n_host)), permute = TRUE)
  # extract dopant points
  rcp_dope = subset(rcp_labeled, marks == "A")
  # extract host points
  #rcp_host = subset(rcp_labeled, marks == "C")

  coords_dope = coords(rcp_dope)
  #win_box = box3(box_2d$xrange, box_2d$yrange, c(0, thickness))
  rcp_box = ppp(x = coords_dope$x, y = coords_dope$y, window = box_2d)

  rcp_G = Gest_nn(rcp_box, correction = "km", k = neighbor_number, r = seq(0, maxGr, length.out = nGr))
  rcp_K =Kest(rcp_box, r = seq(0, maxKr, length.out= nKr))
  summary = list(rrl_G = rcp_G, rrl_K = rcp_K)
  return(summary)
}

#' Perform num_relabs relabelings
#' @export
rcp_simple_relab = function(exp_ppp, num_relabs, thickness, neighbor_number =1, rcp_list, maxGr,
                            maxKr, nGr, nKr, ncores = detectCores()) {

  ## make vectors of maximum radii to take T tests to for K and G
  G_rad_vec = seq(0, maxGr, length.out = (nGr/50) +1)
  G_rad_vec = G_rad_vec[2:length(G_rad_vec)]
  K_rad_vec = seq(0, maxKr, length.out = (nKr/50)+1)
  K_rad_vec = K_rad_vec[2:length(K_rad_vec)]
  nper = num_relabs/floor(length(rcp_list))
  # index for which rcp pattern to use
  ind = rep(1:length(rcp_list), nper)

  ##  calculate envelopes for length(ind) rcp multimer patterns

  cl = makeForkCluster(ncores, outfile = "outfile3")
  rcp_simple = parLapply(cl, 1:length(ind), function(x) {
    func_rcp_simple(exp_ppp, thickness = thickness, rcp_list[[ind[x]]], neighbor_number = neighbor_number,
                    intensity_rcp = 1, maxGr = maxGr, maxKr = maxKr, nGr = nGr, nKr = nKr)
  })


  clusterExport(cl, "rcp_simple", envir = environment())
  # calculate T values for K and G summaries for a variety of max radii
  rcp_simple_T_G = parLapply(cl, G_rad_vec, function(x) {
    Tcalc(rcp_simple, x, func = "G", rmax =maxGr)

  })
  rcp_simple_T_K = parLapply(cl, K_rad_vec, function(x) {
    Tcalc(rcp_simple, x, func = "G", rmax = maxKr)
  })
  # get the 95% CI envelope
  rcp_simple = rrl_averager(rcp_simple, envelope_value = 0.95,  K_cor = "border", G_cor = "km") # use same name to get rid of the older (and enormous) object
  stopCluster(cl)
  return(list(relabs = rcp_simple, T_G = rcp_simple_T_G, T_K = rcp_simple_T_K))
}


### Relabel input point pattern and calcluate summary functions
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
#' @param ... maxKr, nKr, maxGr, nGr, maxGXGHr, maxGXHGr, nGXr
#'
#' @description takes an input binary point pattern, randomly relabels it while maintainig the
#' same ratio of marks, and calculates summary functions
#' @return summary functions listed in `funcs` variable
#' @export
relabeler = function(seed, funcs = c("K", "G", "F", "GXGH"), pattern,
                     dopant_formula, host_formula,
                     K_cor = "trans", G_cor = "km", F_cor = "km",
                     GXGH_cor = "km", GXHG_cor = "km", ...) {
  # set seed to make sure each random relabeling is different
  set.seed(seed)

  #calculate total host and dopant
  dopant_total = sum(marks(pattern)== dopant_formula)
  host_total = sum(marks(pattern)== host_formula)

  # relabel the input pattern, while maintaining the same proportion of host and dopant points
  relabeled = rlabel(pattern,
                         labels = as.factor(c(rep(dopant_formula, dopant_total),
                                              rep(host_formula, host_total))), permute = TRUE)

  # select only the dopant type points
  dopant_relabeled = subset(relabeled, marks == dopant_formula)

  # calculate summary functions
  if (is.pp3(pattern)) {
    if ("K" %in% funcs) {
      if (K_cor == "border") {
        rrl_K = bK3est(dopant_relabeled, rmax = maxKr, nrval =nKr)
      }
      else {
        rrl_K = K3est(dopant_relabeled, rmax = maxKr, nrval = nKr, correction = K_cor)
      }
    }
    if ("G" %in% funcs) {
      rrl_G = G3est(dopant_relabeled, rmax = maxGr, nrval = nGr, correction = G_cor)
    }
    if ("F" %in% funcs) {
      rrl_F = F3est(dopant_relabeled, rmax = maxGr, nrval = nGr, correction = F_cor, vside = vside)
    }
    if ("GXGH" %in% funcs) {
      rrl_GXGH = G3cross(relabeled, i = dopant_formula, j = host_formula,
                         rmax = maxGXGHr, nrval = nGXr, correction = GXGH_cor)
    }
    if ("GXHG" %in% funcs) {
      rrl_GXHG = G3cross(relabeled, i = host_formula, j = dopant_formula,
                         rmax = maxGXHGr, nrval = nGXr, correction =GXHG_cor)
    }
  }
  else if (is.ppp(pattern)) {
    if ("K" %in% funcs) {
      rrl_K = Kest(dopant_relabeled, r = seq(0, maxKr, length.out = nKr), correction = K_cor)
    }
    if ("G" %in% funcs) {
      rrl_G = Gest(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = G_cor)
    }
    if ("F" %in% funcs) {
      rrl_F = Fest(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = F_cor)
    }
    if ("GXGH" %in% funcs) {
      rrl_GXGH = Gcross(relabeled, i = dopant_formula, j = host_formula,
                        r = seq(0, maxGXGHr, length.out = nGXr), correction = GXGH_cor)
    }
    if ("GXHG" %in% funcs) {
      rrl_GXHG = Gcross(relabeled, i = host_formula, j = dopant_formula,
                        r = seq(0, maxGXHGr, length.out = nGXr), correction =GXHG_cor)
    }
  }


  all_funcs = c("K", "G", "F", "GXGH", "GXHG")
  all_funcs %in% funcs
  relabs = c("rrl_K", "rrl_G", "rrl_F", "rrl_GXGH", "rrl_GXHG")
  out =lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) =  relabs[all_funcs %in% funcs]
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
#' of class ppp or pp3.
#' @export
summary_functions = function(funcs = c("K", "G", "F", "GXGH"), pattern, ...,
                             dopant_formula, host_formula,
                             K_cor = "trans", G_cor = "km", F_cor = "km",
                             GXGH_cor = "km", GXHG_cor = "km") {

  # select only the dopant type points
  pattern_dopant = subset(pattern, marks == dopant_formula)

  if (is.pp3(pattern)) {
    # calculate summary functions
    if ("K" %in% funcs) {
      if (K_cor == "border") {
        K = bK3est(pattern_dopant, rmax = maxKr, nrval =nKr)
      }
      else {
        K = K3est(pattern_dopant, rmax = maxKr, nrval = nKr, correction = K_cor)
      }
    }
    if ("G" %in% funcs) {
      G = G3est(pattern_dopant, rmax = maxGr, nrval = nGr, correction = G_cor)
    }
    if ("F" %in% funcs) {
      F = F3est(pattern_dopant, rmax = maxGr, nrval = nGr, correction = F_cor, vside = vside)
    }
    if ("GXGH" %in% funcs) {
      GXGH = G3cross(pattern, i = dopant_formula, j = host_formula,
                     rmax = maxGXGHr, nrval = nGXr, correction = GXGH_cor)
    }
    if ("GXHG" %in% funcs) {
      GXHG = G3cross(pattern, i = host_formula, j = dopant_formula,
                     rmax = maxGXHGr, nrval = nGXr, correction = GXHG_cor)
    }
  }
  else if (is.ppp(pattern)) {
    if ("K" %in% funcs) {
      K = Kest(dopant_relabeled, r = seq(0, maxKr, length.out = nKr), correction = K_cor)
    }
    if ("G" %in% funcs) {
      G = Gest(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = G_cor)
    }
    if ("F" %in% funcs) {
      F = Fest(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = F_cor)
    }
    if ("GXGH" %in% funcs) {
      GXGH = Gcross(relabeled, i = dopant_formula, j = host_formula,
                        r = seq(0, maxGXGHr, length.out = nGXr), correction = GXGH_cor)
    }
    if ("GXHG" %in% funcs) {
      GXHG = Gcross(relabeled, i = host_formula, j = dopant_formula,
                        r = seq(0, maxGXHGr, length.out = nGXr), correction =GXHG_cor)
    }
  }

  all_funcs = c("K", "G", "F", "GXGH", "GXHG")
  relabs = c("K", "G", "F", "GXGH", "GXHG")
  out =lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) =  relabs[all_funcs %in% funcs]
  return(out)
}



#' Plot summary functions
#' @param func Summary function to plot
#' @param observed_values observed summary function value
#' @param relabeled_values theoretical values found by function such as `summary_functions()`
#' @param pattern.colors colors and names for envelope lines
#' @param fill.colors colors and names for envelope fill.  Match names to `pattern.colors`
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description Plot the observed value with envelopes for expected values for summary function
#' @export
summary_plotter = function(func = "K",
                           observed_values, relabeled_values,
                           pattern.colors = c("95.0% AI" = "pink", "Median"  = "black", "Observed" = "blue"),
                           fill.colors = c("95.0% AI" = "pink", "Median"  = "black", "Observed" = "blue"),
                           K_cor = "trans", G_cor = "km", F_cor = "km",
                           GXGH_cor = "km", GXHG_cor = "km",
                           ...)
{
  if (func =="K") {
    baseline = sqrt(relabeled_values$rrl_K$mmean)

    observed = as.data.frame(observed_values$K)
    observed = observed[, K_cor]

    relabeled_values$rrl_K %>% ggplot(aes(x = r,ymin = sqrt(lo) - baseline, ymax = sqrt(hi) - baseline, color = "95.0% AI", fill = "95.0% AI")) +
      geom_ribbon()+
      geom_line(aes( y = sqrt(mmean) - baseline,col = "Median", fill = "Median"), size =2) +
      geom_line(aes(y = sqrt(observed) - baseline, color = "Observed", fill = "Observed"), size = 2) +
      theme(plot.title = element_text(hjust = 0.5, size = 60),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            panel.border = element_rect(linetype = "solid", fill = NA),
            axis.title = element_text(size = 40),
            legend.key.size = unit(20, 'pt'),
            legend.text = element_text(size = 20),
            #legend.title = element_text(size = 20, hjust = 0.5),
            legend.position = c(0.15, 0.80),
            legend.title=element_blank(),
            legend.background = element_rect(fill = "white", color = "black")) +
      guides(color = FALSE) +
      labs(color = "Pattern", fill = "Pattern") + scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab("Radius (nm)") + ylab(expression(tilde(K)[g](r)))

  }

  ## plot for G
  else if (func == "G") {
    baseline = relabeled_values$rrl_G$mmean

    observed = as.data.frame(observed_values$G)
    observed = observed[, G_cor]

    relabeled_values$rrl_G %>% ggplot(aes(x = r, ymin = (lo) - baseline, ymax = (hi) - baseline, color = "95.0% AI", fill = "95.0% AI")) +
      geom_ribbon() +
      geom_line(aes( y = mmean - baseline,col = "Median", fill = "Median"))  +
      geom_line(aes(y = observed - baseline, color = "Observed", fill = "Observed"), size = 2)+
      theme(plot.title = element_text(hjust = 0.5, size = 60),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            panel.border = element_rect(linetype = "solid", fill = NA),
            axis.title = element_text(size = 40),
            legend.key.size = unit(20, 'pt'),
            legend.text = element_text(size = 20),
            #legend.title = element_text(size = 20, hjust = 0.5),
            legend.position = c(0.75, 0.80),
            legend.title=element_blank(),
            legend.background = element_rect(fill = "white", color = "black")) +
      guides(color = FALSE) +
      labs(color = "Pattern", fill = "Pattern") + scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab("Radius (nm)") + ylab(expression(tilde(G)[g](r)))

  }


  ### Plot for F
  else if (func =="F") {
    baseline = relabeled_values$rrl_F$mmean

    observed = as.data.frame(observed_values$F)
    observed = observed[, F_cor]

    relabeled_values$rrl_F%>% ggplot(aes(x = r, ymin = (lo) - baseline, ymax = (hi) - baseline, color = "95.0% AI", fill = "95.0% AI")) +
      geom_ribbon() +
      geom_line(aes( y = (mmean) - baseline,col = "Median", fill = "Median"))  +
      geom_line(aes(y = observed - baseline, color = "Observed", fill = "Observed"), size = 2)+
      theme(plot.title = element_text(hjust = 0.5, size = 60),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            panel.border = element_rect(linetype = "solid", fill = NA),
            axis.title = element_text(size = 40),
            legend.key.size = unit(20, 'pt'),
            legend.text = element_text(size = 20),
            #legend.title = element_text(size = 20, hjust = 0.5),
            legend.position = c(0.75, 0.80),
            legend.title=element_blank(),
            legend.background = element_rect(fill = "white", color = "black")) +
      guides(color = FALSE) +
      labs(color = "Pattern", fill = "Pattern") + scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab("Radius (nm)") + ylab(expression(tilde(F)[g](r)))

  }

  ### Plot for GXGH
  else if (func =="GXGH") {
    baseline = relabeled_values$rrl_GXGH$mmean

    observed = as.data.frame(observed_values$GXGH)
    observed = observed[, GXGH_cor]

    relabeled_values$rrl_GXGH%>% ggplot(aes(x = r, ymin = (lo) - baseline, ymax = (hi) - baseline, color = "95.0% AI", fill = "95.0% AI")) +
      geom_ribbon() +
      geom_line(aes( y = (mmean) - baseline,col = "Median", fill = "Median"))  +
      geom_line(aes(y = observed - baseline, color = "Observed", fill = "Observed"), size = 2)+
      theme(plot.title = element_text(hjust = 0.5, size = 60),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            panel.border = element_rect(linetype = "solid", fill = NA),
            axis.title = element_text(size = 40),
            legend.key.size = unit(20, 'pt'),
            legend.text = element_text(size = 20),
            #legend.title = element_text(size = 20, hjust = 0.5),
            legend.position = c(0.75, 0.80),
            legend.title=element_blank(),
            legend.background = element_rect(fill = "white", color = "black")) +
      guides(color = FALSE) +
      labs(color = "Pattern", fill = "Pattern") + scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab("Radius (nm)") + ylab(expression(tilde(G)[gh](r)))

  }

  ### Plot for GXHG
  else if (func =="GXHG") {
    baseline = relabeled_values$rrl_GXHG$mmean

    observed = as.data.frame(observed_values$GXHG)
    observed = observed[, GXHG_cor]

    relabeled_values$rrl_GXHG%>% ggplot(aes(x = r, ymin = (lo) - baseline, ymax = (hi) - baseline, color = "95.0% AI", fill = "95.0% AI")) +
      geom_ribbon() +
      geom_line(aes( y = (mmean) - baseline,col = "Median", fill = "Median"))  +
      geom_line(aes(y = observed - baseline, color = "Observed", fill = "Observed"), size = 2)+
      theme(plot.title = element_text(hjust = 0.5, size = 60),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            panel.border = element_rect(linetype = "solid", fill = NA),
            axis.title = element_text(size = 40),
            legend.key.size = unit(20, 'pt'),
            legend.text = element_text(size = 20),
            #legend.title = element_text(size = 20, hjust = 0.5),
            legend.position = c(0.75, 0.80),
            legend.title=element_blank(),
            legend.background = element_rect(fill = "white", color = "black")) +
      guides(color = FALSE) +
      labs(color = "Pattern", fill = "Pattern") + scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab("Radius (nm)") + ylab(expression(tilde(G)[hg](r)))

  }
}

#' Tiler function
#' @param pattern input pp3 pattern
#' @param x_tiles,y_tiles,z_tiles number of tiles in each dimension (set equal to 1 for no change)
#' @param z_trim,y_trim,z_trim amount to trim off of each dimension.  Trimming will be applied equall to both sides
#'
#' @description This function takes an input pp3 pattern and tiles it in x, y, and z.  Then, it will trim it
#' according to the `trim` inputs
#' @export
tiler = function(pattern, x_tiles = 1, y_tiles = 1, z_tiles = 1,
                 x_trim = 0, y_trim = 0, z_trim = 0) {
  # get the length of each dimension
  x_length = window(pattern)$domain$xrange[2] - window(pattern)$domain$xrange[1]
  y_length = window(pattern)$domain$yrange[2] - window(pattern)$domain$yrange[1]
  z_length = window(pattern)$domain$zrange[2] - window(pattern)$domain$zrange[1]

  # get new widow
  x_dim =  window(pattern)$domain$xrange[1] + x_length * (x_tiles)
  y_dim =  window(pattern)$domain$yrange[1] + y_length * (y_tiles)
  z_dim =  window(pattern)$domain$zrange[1] + z_length * (z_tiles)
  new_window = as.box3(c(window(pattern)$domain$xrange[1], x_dim),
                       c(window(pattern)$domain$yrange[1], y_dim),
                       c(window(pattern)$domain$zrange[1], z_dim))
  # total number of tiles
  total_tiles = x_tiles * y_tiles * z_tiles
  # total number of points
  total_points = total_tiles * npoints(pattern)

  # save original data
  start_data = pattern$data

  # create matrix will full data
  full_data = pattern$data

  # only add tiles if number of tiles is greater than 1
  if (x_tiles > 1) {
    # new tile for each value x_tiles is above 1
    for (i in 2:x_tiles) {
      # new x_coords are old x_coords with a shift of x_length added
      x_coords =  pattern$data$x + x_length * (i -1)
      dat = start_data
      dat$x = x_coords
      # add the new data to the old data
      full_data = rbind(full_data, dat)
    }
  }
  ## now we tile in y, so we must include the previously tiled data
  start_data = full_data
  if (y_tiles > 1) {
    for (i in 2:y_tiles) {
      y_coords = start_data$y + y_length * (i -1)
      dat = start_data
      dat$y = y_coords
      full_data = rbind(full_data, dat)
    }
  }

  # repeat for z
  start_data = full_data
  if (z_tiles > 1) {
    for (i in 2:z_tiles) {
      z_coords = start_data$z + z_length * (i -1)
      dat = start_data
      dat$z = z_coords
      full_data = rbind(full_data, dat)
    }
  }
  new_pattern = pp3(x =full_data$x, y = full_data$y, z = full_data$z,
                    new_window, marks = full_data$marks)

  ####  get rid of any duplicates that come from tiling with points on the border
  unique_data = as.data.frame.hyperframe(new_pattern$data)
  unique_data = unique(unique_data, margin = 1)
  new_pattern = pp3(x = unique_data$x, y = unique_data$y, z = unique_data$z,
                    new_window, marks = unique_data$marks)
  ##
  ##
  # trim off sides of pattern
  x_min = window(new_pattern)$domain$xrange[1] + x_trim/2
  x_max = window(new_pattern)$domain$xrange[2] - x_trim/2

  y_min = window(new_pattern)$domain$yrange[1] + y_trim/2
  y_max = window(new_pattern)$domain$yrange[2] - y_trim/2

  z_min = window(new_pattern)$domain$zrange[1] + z_trim/2
  z_max = window(new_pattern)$domain$zrange[2] - z_trim/2

  trimmed_window = as.box3(c(x_min, x_max),
                           c(y_min, y_max),
                           c(z_min, z_max))

  trimmed_pattern = pp3(x =new_pattern$data$x, y = new_pattern$data$y,
                        z = new_pattern$data$z,
                        trimmed_window, marks = new_pattern$data$marks)
  trimmed_pattern = subset(trimmed_pattern, subset = trimmed_window)
  return(trimmed_pattern)
}

#' Tile and then relabel
#' @param seed input for `set.seed` function (for reproducibility)
#' @param pp3_full input pp3 pattern for relabeling and tiling
#' @param funcs summary functions to calculate on relabeling
#' @param host_formula,dopant_formula formula for host and dopant marks.  Summary functions
#'  will be calculated from dopant type points
#' @param x_tiles,y_tiles,z_tiles number of tiles in each dimension (set equal to 1 for no change)
#' @param z_trim,y_trim,z_trim amount to trim off of each dimension.  Trimming will be applied equall to both sides
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description This function takes an input pp3 point pattern, relabels it, then
#' tiles the relabeled point pattern.  Then it calculates summary functions on it
tile_relabel = function(seed, pp3_full, funcs = c("K", "G", "F", "GXGH"), ...,
                        host_formula, dopant_formula,
                        x_tiles = 1, y_tiles = 1, z_tiles = 1,
                        x_trim = 0, y_trim = 0, z_trim = 0,
                        K_cor = "border", G_cor = "rs", F_cor = "rs",
                        GXGH_cor = "rs", GXHG_cor = "rs") {

  ## perform relabelings
  # Total number of host and dopant type points
  host_total = sum(pp3_full$data$marks== host_formula)
  dopant_total = sum(pp3_full$data$marks== dopant_formula)
  pp3_dopant = subset(pp3_full, marks == dopant_formula)
  set.seed(seed)

  # relabel the input pattern, while maintaining the same proportion of host and dopant points
  pp3_relabeled = rlabel(pp3_full,
                         labels = as.factor(c(rep(dopant_formula, dopant_total),
                                              rep(host_formula, host_total))), permute = TRUE)

  # select only the dopant type points
  #pp3_dopant_relabeled = subset(pp3_relabeled, marks == dopant_formula)

  pp3_tiled = tiler(pp3_relabeled, x_tiles = x_tiles, y_tiles = y_tiles, z_tiles = z_tiles,
                    x_trim = x_trim, y_trim = y_trim, z_trim = z_trim)

  pp3_dopant_tiled = subset(pp3_tiled, marks == dopant_formula)

  # calculate summary functions
  if ("K" %in% funcs) {
    if (K_cor == "border") {
      rrl_K = bK3est(pp3_dopant_tiled, rmax = maxKr, nrval =nKr)
    }
    else {
      rrl_K = K3est(pp3_dopant_tiled, rmax = maxKr, nrval = nKr, correction = K_cor)
    }
  }
  if ("G" %in% funcs) {
    rrl_G = G3est(pp3_dopant_tiled, rmax = maxGr, nrval = nGr, correction = G_cor)
  }
  if ("F" %in% funcs) {
    rrl_F = F3est(pp3_dopant_tiled, rmax = maxGr, nrval = nGr, correction = F_cor, vside = vside)
  }
  if ("GXGH" %in% funcs) {
    rrl_GXGH = G3cross(pp3_tiled, i = dopant_formula, j = host_formula,
                       rmax = maxGXGHr, nrval = nGXr, correction = GXGH_cor)
  }
  if ("GXHG" %in% funcs) {
    rrl_GXHG = G3cross(pp3_tiled, i = host_formula, j = dopant_formula,
                       rmax = maxGXHGr, nrval = nGXr, correction = GXHG_cor)

  }
  all_funcs = c("K", "G", "F", "GXGH", "GXHG")
  all_funcs %in% funcs
  relabs = c("rrl_K", "rrl_G", "rrl_F", "rrl_GXGH", "rrl_GXHG")
  out =lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) =  relabs[all_funcs %in% funcs]
  return(out)
}






