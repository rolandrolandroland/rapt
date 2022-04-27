#' Rescaler Function
#'
#'#' @description
#' Rescales an input 3D point pattern
#'
#' @details
#' Takes an input 3D point pattern (pp3 object) and scales it to have intensity of new.intensity by multiplying
#' the coordinates of each point as well as the window.
#'
#' @param pattern 3D point pattern (pp3 object) to be scaled
#' @param new.intensity final intensity of point pattern
#'
#' @return pp3.new will preserve all properties of the input \emph{pattern}, except the intensity


rescaler = function(pattern, new.intensity) {
  intensity = spatstat::intensity(pattern)
  factor = (intensity / new.intensity)^(1/3)
  coo = coords(pattern)
  coo.scaled = coo*factor


  win.scaled = box3((domain(pattern)$xrange*factor),
                    (domain(pattern)$yrange*factor),
                    (domain(pattern)$zrange*factor))
  pp3.new = pp3(coo.scaled$x, coo.scaled$y, coo.scaled$z, win.scaled)
}

#' Calculate T value
#'
#' @description
#' Takes input of relabeling object and calculates the T value for either K or G function
#'
#' @details This uses the method from \emph{Baddeley et al.} to calculate the T value for an envelope for either
#' K  (\emph{K3est}) G \emph{G3est} function
#' @references
#' Baddeley, A., Diggle, P. J., Hardegen, A., Lawrence, T., Milne, R. K., & Nair, G. (2014). On tests of spatial pattern based on simulation envelopes. Ecological Monographs, 84(3), 477–489. https://doi.org/10.1890/13-2042.1
Tcalc = function(relabelings, func = "K", rmin = 0, rmax) {
  if (func == "K") {
    mmean =sapply(relabelings, function(x) {
      x$rrl.K$border
    })
    r = relabelings[[1]]$rrl.K$r
    med = apply(mmean, 1, median)
    ind = which(r <=  rmax & r >= rmin)
    interval = rmax / (length(r) -1)
    T = apply(mmean, 2, function(x) {
      T.1= sqrt(x) - sqrt(med)
      T.1 = T.1[ind]
      sum(T.1^2) * interval
    })
    return(T)
  }
  if (func == "G") {
    mmean =sapply(relabelings, function(x) {
      x$rrl.G$km
    })
    r = relabelings[[1]]$rrl.G$r
    med = apply(mmean, 1, median)
    ind = which(r <=  rmax & r >= rmin)
    interval = rmax / (length(r) -1)
    T = apply(mmean, 2, function(x) {
      T.1=x - med
      T.1 = T.1[ind]
      sum(T.1^2) * interval
    })
  }
  return(T)
}

#' 2D relabeling function
relabeler.2d = function(relabelings, envelope.value = .95) {
  mmean =sapply(relabelings, function(x) {
    x$rrl.K$border
  })
  envelope.value =  envelope.value + (1- envelope.value)/2

  ordered = apply(mmean, 1, sort)
  hi.ind = round(length(relabelings) * envelope.value, 0)
  lo.ind = round(length(relabelings) * (1-envelope.value), 0)
  if (lo.ind == 0) {
    lo.ind = 1
  }
  lo = ordered[lo.ind,]
  min = ordered[1,]
  hi = ordered[hi.ind,]
  max = ordered[nrow(ordered),]
  r = relabelings[[1]]$rrl.K$r
  med = apply(mmean, 1, median)
  rrl.K = data.frame(r = r, mmean =med, lo = lo, hi = hi,
                     min = min, max = max,
                     lo.diff = lo - med,
                     hi.diff = hi - med,
                     min.diff = min - med,
                     max.diff = max - med,
                     lo.sqrt.diff = sqrt(lo) - sqrt(med),
                     hi.sqrt.diff = sqrt(hi) - sqrt(med),
                     min.sqrt.diff = sqrt(min) - sqrt(med),
                     max.sqrt.diff = sqrt(max) - sqrt(med))


  # G
  mmean =sapply(relabelings, function(x) {
    x$rrl.G$km
  })
  r = relabelings[[1]]$rrl.G$r
  ordered = apply(mmean, 1, sort)
  lo = ordered[lo.ind,]
  min = ordered[1,]
  hi = ordered[hi.ind,]
  max = ordered[length(relabelings),]
  med = apply(mmean, 1, median)
  rrl.G = data.frame(r = r, mmean =med, lo = lo, hi = hi,
                     min = min, max = max,
                     lo.diff = lo - med,
                     hi.diff = hi - med,
                     min.diff = min - med,
                     max.diff = max - med)

  return(list(rrl.G = rrl.G, rrl.K =rrl.K))
}

#' Simulate 2D Multimers
#' @description Simulate
multimersim = function(exp.ppp,  thickness, group_size = 2,
                       num_neighbors = 6, weights = c(1, 1, 1/4), probs = c(0.6, 0.2, 0.2, 0),
                       rcp.pattern, intensity.rcp = 1, maxGr, maxKr, nGr, nKr) {
  ## make probability vector same length as num_neighbors
  if (num_neighbors > length(probs)) {
    probs = c(probs, rep(0, num_neighbors - length(probs)))
  }

  # get desired intensity
  npoints = exp.ppp$n
  print(npoints)
  box.area = area(box.2d)
  #vol = box.area * thickness
  # intensity.exp = npoints / vol
  # rescale RCP pattern to match physical system (1 point per nm)
  rcp.scaled = rescaler(rcp.pattern, intensity.rcp)

  # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
  rcp.box = subset(rcp.scaled, x > box.2d$xrange[1] & x < box.2d$xrange[2] &
                     y > box.2d$yrange[1] & y < box.2d$yrange[2] &
                     z >0 & z < thickness)
  ## label n/2 points as dopant and the rest as host
  n_irp = npoints
  n_total = length(rcp.box$data$x)
  n_host = n_total - n_irp

  # determine the number of groups of molecules to be present (if multimers, then n_irp/2)
  n_groups = round(n_irp/ group_size,0)
  rcp.labeled = rlabel(rcp.box, labels = c(rep("A",n_groups) , rep("C", length(rcp.box$data$x) - n_groups)), permute = TRUE)
  # extract dopant points
  rcp.dope = subset(rcp.labeled, marks == "A")
  # extract host points
  rcp.host = subset(rcp.labeled, marks == "C")
  # find 6 nearest host to each dopant
  nn.ind= lapply(1:num_neighbors, function(x) {
    nncross(rcp.dope, rcp.host, k = x)
  })

  # get nn distance x, y, and z values for each neighbor
  dist.coords = lapply(1:num_neighbors, function(x) {
    coords(rcp.dope) - coords(rcp.host[nn.ind[[x]][,2]])
  })
  # this just rearranges dist.coords to be grouped by point rather than neighbor order
  holder = c()
  neighbors = lapply(1:nrow(dist.coords[[1]]), function(x) {
    for (i in 1:length(dist.coords)) {
      holder = rbind(holder, dist.coords[[i]][x,])
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
      nn.ind[[which_neighbor[i,j]]][i,2]
    })
  }))

  #duplicate = apply(ind, 2, duplicated)
  #duplicate = duplicated(c(ind))
  duplicate = duplicated(ind, MARGIN = 0)
  duplicate.ind = which(duplicate, arr.ind = TRUE)[,1]
  duplicate.ind = unique(duplicate.ind)
  not.duplicate = ind[!duplicate]
  ## points that are not duplicates
  chosen.points = rcp.host$data[not.duplicate]


  ## host pattern with previously used points removed
  host.unique.pp3 = rcp.host[-ind]
  #unique.points = rcp.host$data[-ind[duplicate],]
  print(paste("first duplicate ", sum(duplicate)))


  ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two dopants )
  # take the 2nd nearest neighbor
  while (sum(duplicate >0) ) {

    # find 6 nearest host to each dopant
    next.nn= lapply(1:num_neighbors, function(x) {
      nncross(rcp.dope, host.unique.pp3, k = x)
    })

    ####
    dist.coords = lapply(1:num_neighbors, function(x) {
      coords(rcp.dope) - coords(host.unique.pp3[next.nn[[x]][,2]])
    })
    ####

    # this just rearranges dist.coords to be grouped by point rather than neighbor order
    holder = c()
    neighbors = lapply(1:nrow(dist.coords[[1]]), function(x) {
      for (i in 1:length(dist.coords)) {
        holder = rbind(holder, dist.coords[[i]][x,])
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
    all.ind= t(sapply(1:nrow(which_neighbor), function(i) {
      sapply(2:(group_size), function(j) {
        next.nn[[which_neighbor[i,j]]][i,2]
      })
    }))

    #ind[duplicate] = all.ind[duplicate]
    #all.ind = c(all.ind)
    # duplicate = c(duplicate)
    # extract just the ones needed to replace duplicates
    next.ind =all.ind[duplicate]

    ## use two duplicates - one that has the index of duplicates in entire all.ind, one that has actual new duplicates
    #new = all.ind[duplicate]
    mat.1 = matrix(FALSE, nrow = nrow(all.ind), ncol = ncol(all.ind))
    mat.1[duplicate] = all.ind[duplicate]
    # zeros = test ==0
    duplicate.mat = duplicated(mat.1, MARGIN = 0)
    duplicate.mat[mat.1==0] = FALSE


    ## get duplicates and their index
    duplicate = duplicated(next.ind, MARGIN = 0)
    duplicate.ind = which(duplicate, arr.ind = TRUE)
    duplicate.ind = unique(duplicate.ind)
    not.duplicate = unique(next.ind)
    ## points that are not duplicates
    new.points = host.unique.pp3$data[not.duplicate]
    chosen.points = rbind(chosen.points, new.points)

    ## host pattern with previously used points removed
    host.unique.pp3 = host.unique.pp3[-next.ind]

    ## change duplicate back to full matrix form
    duplicate = duplicate.mat

    print(paste("loop duplicate ", sum(duplicate)))
  }
  dim(chosen.points)

  # get coordinates of all original points and their selected neighbors
  coords.multimer = rbind(coords(rcp.dope), data.frame(x = chosen.points$x, y = chosen.points$y, z =chosen.points$z))
  dim(coords.multimer)

  # make ppp and caluclate summary functions
  rcp.box = ppp(x = coords.multimer$x, y = coords.multimer$y, window = box.2d)
  rcp.G = Gest(rcp.box, correction = "km", r = seq(0, maxGr, length.out = nGr))
  rcp.K =Kest(rcp.box, r = seq(0, maxKr, length.out= nKr))
  summary = list(rrl.G = rcp.G, rrl.K = rcp.K)
  return(summary)
}

#' Relabel multimers
#'
# multimer
multimers_relab = function (exp.ppp, num.relabs, thickness, rcp.list, maxGr, maxKr,
                            nGr, nKr,
                            group_size = 2, num_neighbors = 6, weights = c(1, 1, 1/4),
                            probs = c(0.6, 0.2, 0.2, 0),
                            ncores = detectCores()) {

  ## make vectors of maximum radii to take T tests to for K and G
  G.rad.vec = seq(0, maxGr, length.out = (nGr/50) +1)
  G.rad.vec = G.rad.vec[2:length(G.rad.vec)]
  K.rad.vec = seq(0, maxKr, length.out = (nKr/50)+1)
  K.rad.vec = K.rad.vec[2:length(K.rad.vec)]

  ## number of relabelings per rcp pattern
  nper = num.relabs/10
  # index for which rcp pattern to use
  ind = rep(1:length(rcp.list), nper)

  cl = makeForkCluster(ncores, outfile = "outfile_test")
  ##  calculate envelopes for length(ind) rcp multimer patterns
  multimers = parLapply(cl, 1:length(ind), function(x) {
    print(paste("start rep  ", x))
    val = multimersim(exp.ppp, thickness = thickness, group_size = group_size,
                      num_neighbors =num_neighbors, weights = weights, probs = probs,
                      rcp.list[[ind[x]]],
                      intensity.rcp = 1, maxGr = maxGr, maxKr = maxKr, nGr = nGr, nKr = nKr)
    print(paste("end rep  ", x))
    val
  })

  clusterExport(cl, "multimers", envir = environment())
  # calculate T values for K and G summaries for a variety of max radii
  multimers.T.G = parLapply(cl, G.rad.vec, function(x) {
    Tcalc(multimers, x, func = "G", rmax =maxGr)
  })
  multimers.T.K = parLapply(cl, K.rad.vec, function(x) {
    Tcalc(multimers, x, func = "K",rmax = maxKr)
  })

  # get the 95% CI envelope
  multimers = relabeler.2d(multimers, envelope.value = 0.95) # use same name to get rid of the older (and enormous) object
  stopCluster(cl)
  return(list(relabs = multimers, T.G = multimers.T.G, T.K = multimers.T.K))
}




## functions for RCP
func.RCP.simple = function(exp.ppp, thickness, rcp.pattern,
                           intensity.rcp = 1, maxGr, maxKr, nGr, nKr) {
  # get desired intensity
  npoints = exp.ppp$n
  print(npoints)
  box.area = area(box.2d)
  #vol = box.area * thickness
  # intensity.exp = npoints / vol
  # rescale RCP pattern to match physical system (1 point per nm)
  rcp.scaled = rescaler(rcp.pattern, intensity.rcp)

  # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
  rcp.box = subset(rcp.scaled, x > box.2d$xrange[1] & x < box.2d$xrange[2] &
                     y > box.2d$yrange[1] & y < box.2d$yrange[2] &
                     z >0 & z < thickness)
  ## label n/2 points as dopant and the rest as host
  n_irp = npoints
  n_total = length(rcp.box$data$x)
  n_host = n_total - n_irp
  rcp.labeled = rlabel(rcp.box, labels = c(rep("A", n_irp) , rep("C", n_host)), permute = TRUE)
  # extract dopant points
  rcp.dope = subset(rcp.labeled, marks == "A")
  # extract host points
  #rcp.host = subset(rcp.labeled, marks == "C")

  coords.dope = coords(rcp.dope)
  #win.box = box3(box.2d$xrange, box.2d$yrange, c(0, thickness))
  rcp.box = ppp(x = coords.dope$x, y = coords.dope$y, window = box.2d)

  rcp.G = Gest(rcp.box, correction = "km", r = seq(0, maxGr, length.out = nGr))
  rcp.K =Kest(rcp.box, r = seq(0, maxKr, length.out= nKr))
  summary = list(rrl.G = rcp.G, rrl.K = rcp.K)
  return(summary)
}

#' Create  RCP Pattern with same intensity and number of points as experimental pattern
func.RCP.simple = function(exp.ppp, thickness, rcp.pattern,
                           intensity.rcp = 1, maxGr, maxKr, nGr, nKr) {
  # get desired intensity
  npoints = exp.ppp$n
  print(npoints)
  box.area = area(box.2d)
  #vol = box.area * thickness
  # intensity.exp = npoints / vol
  # rescale RCP pattern to match physical system (1 point per nm)
  rcp.scaled = rescaler(rcp.pattern, intensity.rcp)

  # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
  rcp.box = subset(rcp.scaled, x > box.2d$xrange[1] & x < box.2d$xrange[2] &
                     y > box.2d$yrange[1] & y < box.2d$yrange[2] &
                     z >0 & z < thickness)
  ## label n/2 points as dopant and the rest as host
  n_irp = npoints
  n_total = length(rcp.box$data$x)
  n_host = n_total - n_irp
  rcp.labeled = rlabel(rcp.box, labels = c(rep("A", n_irp) , rep("C", n_host)), permute = TRUE)
  # extract dopant points
  rcp.dope = subset(rcp.labeled, marks == "A")
  # extract host points
  #rcp.host = subset(rcp.labeled, marks == "C")

  coords.dope = coords(rcp.dope)
  #win.box = box3(box.2d$xrange, box.2d$yrange, c(0, thickness))
  rcp.box = ppp(x = coords.dope$x, y = coords.dope$y, window = box.2d)

  rcp.G = Gest(rcp.box, correction = "km", r = seq(0, maxGr, length.out = nGr))
  rcp.K =Kest(rcp.box, r = seq(0, maxKr, length.out= nKr))
  summary = list(rrl.G = rcp.G, rrl.K = rcp.K)
  return(summary)
}

#' Perform num.relabs relabelings
rcp.simple.relab = function(exp.ppp, num.relabs, thickness, rcp.list, maxGr,
                            maxKr, nGr, nKr, ncores = detectCores()) {

  ## make vectors of maximum radii to take T tests to for K and G
  G.rad.vec = seq(0, maxGr, length.out = (nGr/50) +1)
  G.rad.vec = G.rad.vec[2:length(G.rad.vec)]
  K.rad.vec = seq(0, maxKr, length.out = (nKr/50)+1)
  K.rad.vec = K.rad.vec[2:length(K.rad.vec)]
  nper = num.relabs/10
  # index for which rcp pattern to use
  ind = rep(1:length(rcp.list), nper)

  ##  calculate envelopes for length(ind) rcp multimer patterns

  cl = makeForkCluster(ncores, outfile = "outfile3")
  RCP.simple = parLapply(cl, 1:length(ind), function(x) {
    func.RCP.simple(exp.ppp, thickness = thickness, rcp.list[[ind[x]]],
                    intensity.rcp = 1, maxGr = maxGr, maxKr = maxKr, nGr = nGr, nKr = nKr)
  })


  clusterExport(cl, "RCP.simple", envir = environment())
  # calculate T values for K and G summaries for a variety of max radii
  RCP.simple.T.G = parLapply(cl, G.rad.vec, function(x) {
    Tcalc(RCP.simple, x, func = "G", rmax =maxGr)

  })
  RCP.simple.T.K = parLapply(cl, K.rad.vec, function(x) {
    Tcalc(RCP.simple, x, func = "G", rmax = maxKr)
  })
  # get the 95% CI envelope
  RCP.simple = relabeler.2d(RCP.simple, envelope.value = 0.95) # use same name to get rid of the older (and enormous) object
  stopCluster(cl)
  return(list(relabs = RCP.simple, T.G = RCP.simple.T.G, T.K = RCP.simple.T.K))
}


