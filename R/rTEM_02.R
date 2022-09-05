#' nndist_subset
#'
#' @param pattern  point pattern of type \emph{ppp} or \emph{pp3} to calculate distances from
#' @param pattern2 point pattern of type \emph{ppp} or \emph{pp3} to calculate distances to.
#'  If \emph{NULL} then \emph{pattern} will be used
#' @param window object of class \emph{owin} (for ppp) or \emph{box3} (for pp3).
#' Only points inside the window are used to calculate distances from, but points outside the
#' window are still included as potential neighbors. If \emph{NULL} then it includes the
#' entire domain of \emph{pattern}
#' @param drop_isolated if \emph{TRUE} then points that are closer to a boundary than they
#' are to their \emph{kth} NN will be dropped: their distances are set to \emph{NA}
#' @param k an integer or vector of integers that determines which NN's to calculate
#' @param output outputs a list if "list", outputs a matrix otherwise
#' @description Edge corrected nearest neighbor distribution
#' @details This function calculates the distances to the \emph{kth} nearest neighbors (NN's) for
#' a subset of points, defined as those inside of \emph{window}. All points are still considered
#' when fiding the NN's.
#' @export
nndist_subset = function(pattern, pattern2 = NULL,
                         window = NULL, drop_isolated = FALSE,
                         k =1, output = "list") {

  # if window is not defined
  if (is.null(window) & is.pp3(pattern)) {
    window = domain(pattern)
  }
  if (is.null(window) & is.ppp(pattern)) {
    window = pattern$window
  }

  inside = subset(pattern, inside.boxx(coords(pattern) , w =window))
  if (is.null(pattern2)) {
    pattern2 = pattern
  }

  if (identical(pattern, pattern2)) {
    # first on will just be self - self (0 distance)
    dist = nncross(inside, pattern, k = k+1)
  }
  else {
    dist = nncross(inside, pattern2, k = k)
    ## if pattern is a subset of pattern2
    if (all(dist[,1] == rep(0, length(dist[,1])))) {
      dist = nncross(inside, pattern2, k = k+1)
    }
  }
  if (drop_isolated == TRUE) {
    border_dist = bdist.points(inside)
    head(border_dist)
    ## drop points if closer to border than NN (set to NA)
    out = lapply(1:length(k), function(i) {
      #all distances that are closer to neighbor than border
      distance = sapply(1:length(dist[,i]), function(j){
        if (border_dist[j] >= dist[j,i]) {
          dist[j,i]
        }
        else {
          NA
        }
      })
      # index of each neighbor
      which = sapply(1:length(dist[,i]), function(j){
        if (border_dist[j] >= dist[j,i]) {
          dist[j,i+length(k)]
        }
        else {
          NA
        }
      })
      df = data.frame("distance" = distance, "which" = which)
      names(df) = c(paste("dist", i, sep = "."), paste("which", i, sep = "."))
      df
    })
    names(out) = paste(k, "NN")
  }
  else {
    out = lapply(1:length(k), function(i) {
      distance = dist[,i]
      # index of each neighbor
      which = dist[i+length(k)]
      df = data.frame("distance" = distance, "which" = which)
      names(df) = c(paste("dist", i, sep = "."), paste("which", i, sep = "."))
      df
    })
  }
  if (output == "list") {
    return(out)
  }
    else {
    which_mat = out[[1]][,2]
    dist_mat = out[[1]][,1]
    for (i in 2:length(out)) {
      which_mat = cbind(which_mat, out[[i]][,2])
      dist_mat = cbind(dist_mat, out[[i]][,1])
    }
    out2 = cbind(dist_mat, which_mat)
    colnames(out2)[k] = sapply(k, function(x) {
      paste("dist", x, sep = ".")
    })
    colnames(out2)[(length(k)+1):ncol(out2)] = sapply(k, function(x) {
      paste("which", x, sep = ".")
    })
    return (out2)
  }
}

#' Neighbors in shell
#' @param pattern point pattern of class \emph{ppp} or \emph{pp3}
#' @param type mark defining point type
#' @param window object of class \emph{owin} (for ppp) or \emph{box3} (for pp3).
#' Only points inside the window are used to calculate distances from, but points outside the
#' window are still included as potential neighbors. If \emph{NULL} then it includes the
#' entire domain of \emph{pattern}
#' @param k_vec an integer or vector of integers that determines which NN's to calculate
#' @param drop_isolated if \emph{TRUE} then points that are closer to a boundary than they
#' are to their \emph{kth} NN will be dropped: their distances are set to \emph{NA}
#' @description Compare distribution of marks to binomial distribution
#' @details For each k in \emph{k_vec}, calculate how many of the points in
#' \emph{pattern} with mark \emph{type} have  0, 1, .. \emph{k} points of mark \emph{type} in their \emph{k}
#' nearest neighbors.
#' This can then be compared to the expected values of a binomial distribution with n = k and
#' p = fraction of points in \emph{pattern} that have mark \emph{type}
#' @export
neighbors = function(pattern, type = NULL,
                     window = NULL, k_vec = 1,
                     drop_isolated = TRUE,
                     ...) {
  ## k_vec must have every integer for this to work
  k_vec = 1:max(k_vec)

  # subset pattern to those with marks type
  pattern_subset = subset(pattern, marks %in% type)

  # for each k in k_vec, get the ditance from each point to its kth nearest neighbor and
  # and the identity of that NN
    # NA means that the point is closer to border than that NN (only used if drop_isolated is TRUE)
  nndist_g_all = nndist_subset(pattern = pattern_subset, pattern2 = pattern,
                               window = window,
                               drop_isolated = drop_isolated, k = k_vec)

  # now get the identity of each of those nearest neighbors
  g_all_marks = sapply(nndist_g_all, function(i) {
    marks(pattern)[i$which]
  })
  g_all_marks = g_all_marks[complete.cases(g_all_marks),]
  ## for each point, this returns the number of points within that shell that are same type
  num_neighbors = as.data.frame(sapply(k_vec, function(shell_size) {
    apply(g_all_marks, 1, function(i) {
     sum(i[1:shell_size] %in% type)
    })
  }))
  shell_long = unlist(lapply(k_vec, rep, nrow(num_neighbors)))
  num = c()
  for (j in 1:ncol(num_neighbors)) {
    num = c(num, num_neighbors[,j])
  }
  num_neighbors_long = data.frame(num = num,
                                  shell= shell_long)
  sum = sapply(k_vec, function(k) {
    vec = num_neighbors_long %>% filter(shell ==k)
    #print(vec)
    sapply(0:length(k_vec), function(inner) {
      sum(vec$num == (inner))
    })
  })
  sum = as.data.frame(sum)
  rownames(sum)= c(0:length(k_vec))
  colnames(sum) = k_vec
  sum
}

#' nn_binom
#' @param k_vec vector of integers to use as \emph{n} values for binomial distribution
#' @param p probability of success. Should equal the concentration of the point pattern that
#' this will be compared to using the \emph{neighbors()} function
#' @description calculates binomial distributions with \emph{n = k}
#' for each \emph{k} in \emph{k_vec} and probability of success \emph{p}
#' @export
nn_binom = function(k_vec, p) {
  sapply(k_vec, function(shell) {
    dbinom(0:max(k_vec), shell, pcp)
  })
}

#' Multimer Simulation
#' @param dopant_pattern point pattern of class \emph{ppp} or \emph{pp3}.  The final multimer
#' pattern will match this pattern in class, intensity, and domain.
#' @param upp point pattern of class \emph{ppp} or \emph{pp3} to use as underlying point pattern.
#' Multimers will be selected from these points.
#' @param min_thick if \emph{dopant_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the smallest z value to keep before collapsing into 2d.
#' @param max_thick if \emph{dopant_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the largest z value to keep before collapsing into 2d.
#' @param group_size size of clusters.  If equal to 1, then all points will be independent
#'  of each other
#' @param num_neighbors number of nearest neighbors to select from when
#' forming dimers, trimers, etc..
#' @param weights vector of length equal to number of dimensions in \emph{upp}. the weighted distance to each of \emph{num_neighbors} nearest neighbors
#' is calculated using \eqn{\sqrt{w_1 x^2 + w_2 y^2 + w_3 z^2}},
#'  where \emph{weights} = (\eqn{w_1, w_2, w_3})
#' @param probs vector of probabilities.  Should sum to \emph{group_size}-1.
#' For \eqn{probs = c(p_1, p_2, p_3, p_4)}, the probability of the first NN being
#' selected in \eqn{p_1}, the probability of the second is \eqn{p_2}, and so on
#' @param intensity_upp the \emph{upp} will be rescaled to have this intensity before the
#' marks are assigned.
#' @description Simulates multimers (groups of two/dimers, three/trimers, etc.) in
#' underlying point pattern (upp) to match domain and number of points in \emph{dopant_pattern}.
#' @details algorithm steps:
#'  \itemize{
#'  \item{"Step 1"}{rescale \emph{upp} to match \emph{intensity_upp}}
#'  \item{"Step 2"}{Select points in the scaled \emph{upp} that are
#'  inside the domain of \emph{dopant_pattern}. }
#' }
#'
#' @export
multimersim = function(dopant_pattern, upp, min_thick = 0, max_thick = 10, group_size = 2,
                            num_neighbors = 6, weights = c(1, 1, 1), probs = c(1, 0, 0, 0),
                            intensity_upp = NULL) {
  ## make probability vector same length as num_neighbors
  if (num_neighbors > length(probs)) {
    probs = c(probs, rep(0, num_neighbors - length(probs)))
  }
  if (is.null(intensity_upp)) {
    intensity_upp = sum(spatstat.geom::intensity(upp))
  }
  # get desired intensity
  npoints = npoints(dopant_pattern)
  print(npoints)

  # rescale UPP pattern to match physical system (1 point per nm)
  upp_scaled = rTEM::rescaler(upp, intensity_upp)

  # case 1 and 2:  dopant_pattern is 2d
  if (is.ppp(dopant_pattern)) {
    box_2d = dopant_pattern$window
    box_area = spatstat.geom::area(box_2d)

    # case 1 UPP pattern is 3d
    if (is.pp3(upp)) {
      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2] &
                         z >=min_thick & z <= max_thick)
      ## label n/2 points as dopant and the rest as host
      n_irp = npoints
      n_total = npoints(upp_box)
      n_host = n_total - n_irp

      # determine the number of groups of molecules to be present (if multimers, then n_irp/2)
      n_groups = round(n_irp/ group_size,0)
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
      # extract dopant points
      upp_dope = subset(upp_labeled, marks == "A")
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")
      if (group_size == 1) {
        multimer_box = ppp(x = upp_dope$data$x, y = upp_dope$data$y, window = box_2d)
        return(multimer_box)
        next
      }
      # find 6 nearest host to each dopant
      nn_ind= lapply(1:num_neighbors, function(x) {
        nncross(upp_dope, upp_host, k = x)
      })
      # get nn distance x, y, and z values for each neighbor
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_dope) - coords(upp_host[nn_ind[[x]][,2]])
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
          sqrt(weights[1]*neighbors[[outer]]$x[inner]^2 + weights[2]*neighbors[[outer]]$y[inner]^2 + weights[3]*neighbors[[outer]]$z[inner]^2)
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
      chosen_points = coords(upp_host)[not_duplicate,]


      ## host pattern with previously used points removed
      host_unique = upp_host[-ind]
      #unique_points = upp_host$data[-ind[duplicate],]
      print(paste("first duplicate ", sum(duplicate)))


      ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two dopants )
      # take the 2nd nearest neighbor
      while (sum(duplicate >0) ) {

        # find 6 nearest host to each dopant
        next_nn= lapply(1:num_neighbors, function(x) {
          nncross(upp_dope, host_unique, k = x)
        })

        ####
        dist_coords = lapply(1:num_neighbors, function(x) {
          coords(upp_dope) - coords(host_unique[next_nn[[x]][,2]])
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

        # rank each of the nearest neighbors based upon distance and weights (distance will be ONLY x-y distance)
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


        # extract just the ones needed to replace duplicates
        next_ind =all_ind[duplicate]

        ## use two duplicates - one that has the index of duplicates in entire all_ind, one that has actual new duplicates
        mat_1 = matrix(FALSE, nrow = nrow(all_ind), ncol = ncol(all_ind))
        mat_1[duplicate] = all_ind[duplicate]
        duplicate_mat = duplicated(mat_1, MARGIN = 0)
        duplicate_mat[mat_1==0] = FALSE


        ## get duplicates and their index
        duplicate = duplicated(next_ind, MARGIN = 0)
        duplicate_ind = which(duplicate, arr.ind = TRUE)
        duplicate_ind = unique(duplicate_ind)
        not_duplicate = unique(next_ind)
        ## points that are not duplicates
        new_points = coords(host_unique)[not_duplicate,]
        chosen_points = rbind(chosen_points, new_points)

        ## host pattern with previously used points removed
        host_unique = host_unique[-next_ind]

        ## change duplicate back to full matrix form
        duplicate = duplicate_mat

        print(paste("loop duplicate ", sum(duplicate)))
      }
      dim(chosen_points)

      # get coordinates of all original points and their selected neighbors
      coords_multimer = rbind(coords(upp_dope)[,c(1,2)], data.frame(x = chosen_points$x, y = chosen_points$y))
      dim(coords_multimer)

      # make ppp and caluclate summary functions
      multimer_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
    }
    ### END CASE 1


    # case 2 UPP is 2d
    else if (is.ppp(upp)) {
      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2])
      ## label n/2 points as dopant and the rest as host
      n_irp = npoints
      n_total = npoints(upp_box)
      n_host = n_total - n_irp

      # determine the number of groups of molecules to be present (if multimers, then n_irp/2)
      n_groups = round(n_irp/ group_size,0)
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
      # extract dopant points
      upp_dope = subset(upp_labeled, marks == "A")
      if (group_size == 1) {
        multimer_box = ppp(x = upp_dope$x, y = upp_dope$y, window = box_2d)
        return(multimer_box)
        next
      }
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")
      # find 6 nearest host to each dopant
      nn_ind= lapply(1:num_neighbors, function(x) {
        nncross(upp_dope, upp_host, k = x)
      })
      # get nn distance x, y, and z values for each neighbor
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_dope) - coords(upp_host[nn_ind[[x]][,2]])
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
          sqrt(weights[1]*neighbors[[outer]]$x[inner]^2 + weights[2]*neighbors[[outer]]$y[inner]^2)
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
      chosen_points = coords(upp_host)[not_duplicate,]


      ## host pattern with previously used points removed
      host_unique = upp_host[-ind]
      #unique_points = upp_host$data[-ind[duplicate],]
      print(paste("first duplicate ", sum(duplicate)))


      ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two dopants )
      # take the 2nd nearest neighbor
      while (sum(duplicate >0) ) {

        # find 6 nearest host to each dopant
        next_nn= lapply(1:num_neighbors, function(x) {
          nncross(upp_dope, host_unique, k = x)
        })

        ####
        dist_coords = lapply(1:num_neighbors, function(x) {
          coords(upp_dope) - coords(host_unique[next_nn[[x]][,2]])
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

        # rank each of the nearest neighbors based upon distance and weights (distance will be ONLY x-y distance)
        # rank will be from smallest to largest
        ranks = lapply(1:length(neighbors), function(outer) {
          distances = sapply(1:nrow(neighbors[[outer]]), function(inner) {
            sqrt(weights[1]*neighbors[[outer]]$x[inner]^2 + weights[2]*neighbors[[outer]]$y[inner]^2)
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


        # extract just the ones needed to replace duplicates
        next_ind =all_ind[duplicate]

        ## use two duplicates - one that has the index of duplicates in entire all_ind, one that has actual new duplicates
        mat_1 = matrix(FALSE, nrow = nrow(all_ind), ncol = ncol(all_ind))
        mat_1[duplicate] = all_ind[duplicate]
        duplicate_mat = duplicated(mat_1, MARGIN = 0)
        duplicate_mat[mat_1==0] = FALSE


        ## get duplicates and their index
        duplicate = duplicated(next_ind, MARGIN = 0)
        duplicate_ind = which(duplicate, arr.ind = TRUE)
        duplicate_ind = unique(duplicate_ind)
        not_duplicate = unique(next_ind)
        ## points that are not duplicates
        new_points = coords(host_unique)[not_duplicate,]
        chosen_points = rbind(chosen_points, new_points)

        ## host pattern with previously used points removed
        host_unique = host_unique[-next_ind]

        ## change duplicate back to full matrix form
        duplicate = duplicate_mat

        print(paste("loop duplicate ", sum(duplicate)))
      }
      dim(chosen_points)

      # get coordinates of all original points and their selected neighbors
      coords_multimer = rbind(coords(upp_dope), data.frame(x = chosen_points$x, y = chosen_points$y))
      dim(coords_multimer)

      # make ppp and caluclate summary functions
      multimer_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
    }
    return(multimer_box)
    ## END CASE 2
  }

  # case 3: dopant_pattern is 3d and upp is 3d
  else if (is.pp3(dopant_pattern)) {
    box_3d = domain(dopant_pattern)
    box_area = spatstat.geom::volume(box_3d)
    # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
    upp_box = subset(upp_scaled, x >= box_3d$xrange[1] & x <= box_3d$xrange[2] &
                       y >= box_3d$yrange[1] & y <= box_3d$yrange[2] &
                       z >= box_3d$zrange[1] & z <= box_3d$zrange[2])
    ## label n/2 points as dopant and the rest as host
    n_irp = npoints
    n_total = npoints(upp_box)
    n_host = n_total - n_irp

    # determine the number of groups of molecules to be present (if multimers, then n_irp/2)
    n_groups = round(n_irp/ group_size,0)
    upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
    # extract dopant points
    upp_dope = subset(upp_labeled, marks == "A")

    if (group_size == 1) {
      multimer_box = pp3(x = upp_dope$data$x, y = upp_dope$data$y,z = upp_dope$data$z, window = box_3d)
      ind = match(do.call("paste", coords(multimer_box)), do.call("paste", coords(upp_labeled)))
      marks(upp_labeled) = "H"
      marks(upp_labeled)[ind] = "G"
      marks(upp_labeled) = as.factor(marks(upp_labeled))
      return(upp_labeled)
      next
    }
    # extract host points
    upp_host = subset(upp_labeled, marks == "C")
    # find 6 nearest host to each dopant
    nn_ind= lapply(1:num_neighbors, function(x) {
      nncross(upp_dope, upp_host, k = x)
    })
    # get nn distance x, y, and z values for each neighbor
    dist_coords = lapply(1:num_neighbors, function(x) {
      coords(upp_dope) - coords(upp_host[nn_ind[[x]][,2]])
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
        sqrt(weights[1]*neighbors[[outer]]$x[inner]^2 + weights[2]*neighbors[[outer]]$y[inner]^2 + weights[3]*neighbors[[outer]]$z[inner]^2)
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
    chosen_points = coords(upp_host)[not_duplicate,]


    ## host pattern with previously used points removed
    host_unique = upp_host[-ind]
    #unique_points = upp_host$data[-ind[duplicate],]
    print(paste("first duplicate ", sum(duplicate)))


    ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two dopants )
    # take the 2nd nearest neighbor
    while (sum(duplicate >0) ) {

      # find 6 nearest host to each dopant
      next_nn= lapply(1:num_neighbors, function(x) {
        nncross(upp_dope, host_unique, k = x)
      })

      ####
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_dope) - coords(host_unique[next_nn[[x]][,2]])
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

      # rank each of the nearest neighbors based upon distance and weights (distance will be ONLY x-y distance)
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


      # extract just the ones needed to replace duplicates
      next_ind =all_ind[duplicate]

      ## use two duplicates - one that has the index of duplicates in entire all_ind, one that has actual new duplicates
      mat_1 = matrix(FALSE, nrow = nrow(all_ind), ncol = ncol(all_ind))
      mat_1[duplicate] = all_ind[duplicate]
      duplicate_mat = duplicated(mat_1, MARGIN = 0)
      duplicate_mat[mat_1==0] = FALSE


      ## get duplicates and their index
      duplicate = duplicated(next_ind, MARGIN = 0)
      duplicate_ind = which(duplicate, arr.ind = TRUE)
      duplicate_ind = unique(duplicate_ind)
      not_duplicate = unique(next_ind)
      ## points that are not duplicates
      new_points = coords(host_unique)[not_duplicate,]
      chosen_points = rbind(chosen_points, new_points)

      ## host pattern with previously used points removed
      host_unique = host_unique[-next_ind]

      ## change duplicate back to full matrix form
      duplicate = duplicate_mat

      print(paste("loop duplicate ", sum(duplicate)))
    }
    dim(chosen_points)

    # get coordinates of all original points and their selected neighbors
    coords_multimer = rbind(coords(upp_dope),
                            data.frame(x = chosen_points$x, y = chosen_points$y, z = chosen_points$z))
    dim(coords_multimer)

    # make ppp and caluclate summary functions
    multimer_box = pp3(x = coords_multimer$x, y = coords_multimer$y, z = coords_multimer$z, window = box_3d)
    ind = match(do.call("paste", coords(multimer_box)), do.call("paste", coords(upp_box)))
    marks(upp_box) = "H"
    marks(upp_box)[ind] = "G"
    marks(upp_box) = as.factor(marks(upp_box))
    return(upp_box)
  }
}

#' pp3 downsizer
#' @export
pp3_downsizer = function(ranged_pos, downsize, side = "both") {
  if (side == "both") { ## remove from both sides
    win_new = box3(c(min(ranged_pos$x) + downsize, max(ranged_pos$x) - downsize),
                   c(min(ranged_pos$y) +downsize, max(ranged_pos$y)- downsize),
                   c(min(ranged_pos$z) + downsize, max(ranged_pos$z)- downsize))
    # one with all points
    ranged_pos_new = ranged_pos[ranged_pos$x < max(ranged_pos$x) - downsize &  ranged_pos$x >min(ranged_pos$x) + downsize &
                                  ranged_pos$y < max(ranged_pos$y) - downsize &  ranged_pos$y >min(ranged_pos$y) + downsize &
                                  ranged_pos$z < max(ranged_pos$z) - downsize &  ranged_pos$z >min(ranged_pos$z) + downsize,]
  }
  else if (side == "negative") { # just remove from min
    win_new = box3(c(min(ranged_pos$x) + downsize, max(ranged_pos$x)),
                   c(min(ranged_pos$y) +downsize, max(ranged_pos$y)),
                   c(min(ranged_pos$z) + downsize, max(ranged_pos$z)))
    # one with all points
    ranged_pos_new = ranged_pos[ranged_pos$x < max(ranged_pos$x) &  ranged_pos$x >min(ranged_pos$x) + downsize &
                                  ranged_pos$y < max(ranged_pos$y)  &  ranged_pos$y >min(ranged_pos$y) + downsize &
                                  ranged_pos$z < max(ranged_pos$z)  &  ranged_pos$z >min(ranged_pos$z) + downsize,]
  }
  else if (side == "positive") { # just remove from max
    win_new = box3(c(min(ranged_pos$x) , max(ranged_pos$x) - downsize),
                   c(min(ranged_pos$y), max(ranged_pos$y)- downsize),
                   c(min(ranged_pos$z), max(ranged_pos$z)- downsize))
    # one with all points
    ranged_pos_new = ranged_pos[ranged_pos$x < max(ranged_pos$x) - downsize &  ranged_pos$x >min(ranged_pos$x)  &
                                  ranged_pos$y < max(ranged_pos$y) - downsize &  ranged_pos$y >min(ranged_pos$y)  &
                                  ranged_pos$z < max(ranged_pos$z) - downsize &  ranged_pos$z >min(ranged_pos$z) ,]
  }
  else {
    print("side variable must be either 'both', 'positive', or 'negative'")
    return(NULL)
  }

  pp3_new = createSpat(ranged_pos_new, win = win_new)
  marks(pp3_new) = as.factor(ranged_pos_new$mark)
  return(pp3_new)
}

#' Collpase pp3 into ppp
#' @description function that collapses a pp3 into a ppp
#' @export
pp3_collapse = function(input_data, dim = "z", output = "ppp") {
  all_dims = c("x", "y", "z")
  output_dims = all_dims[!(all_dims %in% dim)]
  if (first(class(input_data)) == "pp3")
  {
    data = input_data$data
    mark_data = marks(input_data)
    # get all ranges
    xrange = input_data$domain$xrange
    yrange = input_data$domain$yrange
    zrange = input_data$domain$zrange

    ranges = c("xrange", "yrange", "zrange")

    # select the ones that we are keeping
    output_range = ranges[!(all_dims %in% dim)]
    range1 = lapply(output_range, get, envir = environment())[[1]]
    range2 = lapply(output_range, get, envir = environment())[[2]]


    win = owin(range1, range2)
  }
  ## otherwise it should be dataframe
  else
  {
    data = input_data
    mark_data = input_data$mark

    xrange = range(data$x)
    yrange = range(data$y)
    zrange = range(data$z)
    ranges = c("xrange", "yrange", "zrange")

    # select the ones that we are keeping
    output_range = ranges[!(all_dims %in% dim)]
    range1 = lapply(output_range, get, envir = environment())[[1]]
    range2 = lapply(output_range, get, envir = environment())[[2]]
    win = owin(range1, range2)

  }

  xvals = data$x
  yvals = data$y
  zvals = data$z
  all_vals = c("xvals", "yvals", "zvals")
  output_vals = all_vals[!(all_dims %in% dim)]
  vals_1 = lapply(output_vals, get, envir = environment())[[1]]
  vals_2 = lapply(output_vals, get, envir = environment())[[2]]

  if (output == "ppp")
  {

    new_data = ppp(vals_1, vals_2, window = win)
    marks(new_data) = mark_data

  }
  else
  {
    new_data = data.frame(vals_1, vals_2)
    colnames(new_data) = output_dims
    new_data$mark = mark_data

  }
  return(new_data)
}
