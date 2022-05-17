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
#' @description
#' Takes input of relabeling object and calculates the T value for either K or G function
#'
#' @details This uses the method from \emph{Baddeley et al.} to calculate the T value for an envelope for either
#' K  (\emph{K3est}) G \emph{G3est} function
#' @references
#' Baddeley, A., Diggle, P. J., Hardegen, A., Lawrence, T_, Milne, R. K., & Nair, G. (2014). On tests of spatial pattern based on simulation envelopes. Ecological Monographs, 84(3), 477–489. https://doi.org/10.1890/13-2042.1
Tcalc = function(relabelings, func = "K", rmin = 0, rmax) {
  if (func == "K") {
    mmean =sapply(relabelings, function(x) {
      x$rrl_K$border
    })
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
    mmean =sapply(relabelings, function(x) {
      x$rrl_G$km
    })
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

#' 2D relabeling function
relabeler_2d = function(relabelings, envelope_value = .95) {
  mmean =sapply(relabelings, function(x) {
    x$rrl_K$border
  })
  envelope_value =  envelope_value + (1- envelope_value)/2

  ordered = apply(mmean, 1, sort)
  hi_ind = round(length(relabelings) * envelope_value, 0)
  lo_ind = round(length(relabelings) * (1-envelope_value), 0)
  if (lo_ind == 0) {
    lo_ind = 1
  }
  lo = ordered[lo_ind,]
  min = ordered[1,]
  hi = ordered[hi_ind,]
  max = ordered[nrow(ordered),]
  r = relabelings[[1]]$rrl_K$r
  med = apply(mmean, 1, median)
  rrl_K = data.frame(r = r, mmean =med, lo = lo, hi = hi,
                     min = min, max = max,
                     lo_diff = lo - med,
                     hi_diff = hi - med,
                     min_diff = min - med,
                     max_diff = max - med,
                     lo_sqrt_diff = sqrt(lo) - sqrt(med),
                     hi_sqrt_diff = sqrt(hi) - sqrt(med),
                     min_sqrt_diff = sqrt(min) - sqrt(med),
                     max_sqrt_diff = sqrt(max) - sqrt(med))


  # G
  mmean =sapply(relabelings, function(x) {
    x$rrl_G$km
  })
  r = relabelings[[1]]$rrl_G$r
  ordered = apply(mmean, 1, sort)
  lo = ordered[lo_ind,]
  min = ordered[1,]
  hi = ordered[hi_ind,]
  max = ordered[length(relabelings),]
  med = apply(mmean, 1, median)
  rrl_G = data.frame(r = r, mmean =med, lo = lo, hi = hi,
                     min = min, max = max,
                     lo_diff = lo - med,
                     hi_diff = hi - med,
                     min_diff = min - med,
                     max_diff = max - med)

  return(list(rrl_G = rrl_G, rrl_K =rrl_K))
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
# multimer
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
  multimers = relabeler_2d(multimers, envelope_value = 0.95) # use same name to get rid of the older (and enormous) object
  stopCluster(cl)
  return(list(relabs = multimers, T_G = multimers_T_G, T_K = multimers_T_K))
}




## functions for RCP
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
  rcp_simple = relabeler_2d(rcp_simple, envelope_value = 0.95) # use same name to get rid of the older (and enormous) object
  stopCluster(cl)
  return(list(relabs = rcp_simple, T_G = rcp_simple_T_G, T_K = rcp_simple_T_K))
}


