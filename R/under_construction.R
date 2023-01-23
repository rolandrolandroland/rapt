#' Multimer Simulation development version
#' @param guest_pattern point pattern of class \emph{ppp} or \emph{pp3}.  The final multtimer
#' pattern will match this pattern in class, intensity, and domain.  If this is left as NULL, then
#' the domain will match that of \emph{upp}, will be of the class specified in \emph{output},
#' and have number of guests specified in \emph{n_guest}
#' @param upp point pattern of class \emph{ppp} or \emph{pp3} to use as underlying point pattern.
#' Multimers will be selected from these points.
#' @param output leave as \emph{"guest pattern type"} if specifying \emph{guest_pattern}.  Otherwise, set to \emph{"ppp"} for
#' 2D output or \emph{pp3} for 3D output
#' @param n_guests leave as \emph{NA} if specifying \emph{UPP}. Otherwise, set to
#' integer for number of guests in multimer pattern
#' @param min_thick if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the smallest z value to keep before collapsing into 2d.
#' @param max_thick if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the largest z value to keep before collapsing into 2d.
#' @param ztrim a numeric.  This will trim this amount from both top and bottom (positive and negative z)
#' after multimers are generated and before pp3 pattern is collapsed into ppp.
#' Only applies if \emph{upp} is 3D (\emph{pp3})
#' @param group_size size of clusters.  If equal to 1, then all points will be independent
#'  of each other
#' @param num_neighbors number of nearest neighbors to select from when
#' forming dimers, trimers, etc..
#' @param weights vector of length equal to number of dimensions in \emph{upp}. the weighted distance to each of \emph{num_neighbors} nearest neighbors
#' is calculated using \eqn{\sqrt{w_1 x^2 + w_2 y^2 + w_3 z^2}},
#'  where \emph{weights} = (\eqn{w_1, w_2, w_3}). Set to \emph{c(1, 1, 0)} for vertical dimers.
#' @param probs vector of probabilities.  Should sum to \emph{group_size}-1.
#' For \eqn{probs = c(p_1, p_2, p_3, p_4)}, the probability of the first NN being
#' selected in \eqn{p_1}, the probability of the second is \eqn{p_2}, and so on
#' @param intensity_upp the \emph{upp} will be rescaled to have this intensity before the
#' marks are assigned. Leave as \emph{NA} to use \emph{upp} as it is
#' @description Under construction. See \code{\link{multimersim}} for stable version. Simulates multimers (groups of two/dimers, three/trimers, etc.) in
#' underlying point pattern (upp) to match domain and number of points in \emph{guest_pattern}.
#' @details algorithm steps:
#'  \itemize{
#'  \item{Step 1:} {rescale \emph{upp} to match \emph{intensity_upp}}
#'  \item{Step 2:} {Select points in the scaled \emph{upp} that are
#'  inside the domain of \emph{guest_pattern}. }
#'  \item{Step 3:} {Determine number of guest groups or clusters (for dimers, this is number of guests / 2)
#'  and assign this many points in the scaled subsetted UPP to be guests.
#'   These are the "centroids"}
#'  \item{Step 4:} {Take the \emph{num_neighbors} closest point to each guest}
#'  \item{Step 5:} {Rank each neighbor by weighted distance to centroid}
#'  \item{Step 6:} {Using the probabilities in \emph{probs}, select which neighbors
#'   are to be guests (so that the cluster size is now equal to \emph{group_size})}
#'  \item{Step 7:} {For any duplicates, redo process so that correct number of guests are present}
#'  \item{Step 8:} {If \emph{guest_pattern} is 2D and \emph{UPP} is 3D, remove Z coordinate to
#'  make new pattern 2D}
#'}
#'
#' @export

multimersim_dev = function(guest_pattern = NULL, upp, output = "guest pattern type", n_guests = NA,
                           min_thick = NA, max_thick = NA, ztrim = 0,
                           group_size = 2, num_neighbors = 6,
                           weights = c(1, 1, 1), probs = c(1, 0, 0, 0),
                           intensity_upp = NA) {

  ## check input parameters
  if(is.null(guest_pattern)) {
    if (output == "guest pattern type" || is.na(n_guests)) {
      stop("guest_pattern is not defined: output must be either `ppp` or `pp3` and
         n_guests must be a numeric")
    }
    ## check that guest has fewer points than UPP
    if (n_guests >= spatstat.geom::npoints(upp)) {
      stop("n_guests must be smaller than number of points in upp")
    }
  }
  else if (spatstat.geom::npoints(guest_pattern) >= spatstat.geom::npoints(upp)) {
    stop("guest pattern must have fewer points than upp")
  }

  ## make probability vector same length as num_neighbors
  if (num_neighbors > length(probs)) {
    probs = c(probs, rep(0, num_neighbors - length(probs)))
  }
  if (is.na(intensity_upp)) {
    intensity_upp = sum(spatstat.geom::intensity(upp))
  }

  # rescale UPP pattern to match physical system (1 point per nm)
  upp_scaled = rTEM::rescaler(upp, intensity_upp)

  # case 1 and 2:  guest_pattern is 2d
  if (spatstat.geom::is.ppp(guest_pattern) || output == "ppp") {

    if (is.null(guest_pattern) && spatstat.geom::is.pp3(upp)) {
      box_2d = as.owin(c(upp$domain$xrange,
                         upp$domain$yrange))

    }
    else if (is.null(guest_pattern) && spatstat.geom::is.ppp(upp)) {
      box_2d = upp$window
    }

    else {
      box_2d = guest_pattern$window
      n_guests = npoints(guest_pattern)
    }

    box_area = spatstat.geom::area(box_2d)

    # case 1 UPP pattern is 3d
    if (spatstat.geom::is.pp3(upp)) {

      if (is.na(max_thick)) {
        max_thick = max(upp$domain$zrange)
      }
      if (is.na(min_thick)) {
        min_thick = min(upp$domain$zrange)
      }
      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2] &
                         z >=min_thick - ztrim & z <= max_thick + ztrim)

      ## If using ztrim, then must adjust so that trimmed box has correct number of points
      # if not using ztrim, it will be zero and full/final will just be 1
      full_volume = abs(box_2d$xrange[2] - box_2d$xrange[1]) *
        abs(box_2d$yrange[2] - box_2d$yrange[1]) *
        (max_thick - min_thick + ztrim *2)
      final_volume = abs(box_2d$xrange[2] - box_2d$xrange[1]) *
        abs(box_2d$yrange[2] - box_2d$yrange[1]) *
        (max_thick - min_thick)

      # npoints will be scaled by volume ratios
      ## label n/2 points as guest and the rest as host
      n_guest = n_guests * full_volume/final_volume
      n_total = npoints(upp_box) #* full_volume/final_volume
      n_host = n_total - n_guest

      # determine the number of groups of molecules to be present (if dimers, then n_guest/2)
      n_groups = round(n_guest/ group_size,0)
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
      # extract guest points
      upp_guest = subset(upp_labeled, marks == "A")
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")
      if (group_size == 1) {
        multimer_box = ppp(x = upp_guest$data$x, y = upp_guest$data$y, window = box_2d)
        return(multimer_box)
        next
      }
      # find 6 nearest host to each guest
      nn_ind= lapply(1:num_neighbors, function(x) {
        nncross(upp_guest, upp_host, k = x)
      })
      # get nn distance x, y, and z values for each neighbor
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_guest) - coords(upp_host[nn_ind[[x]][,2]])
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


      ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two guests )
      # take the 2nd nearest neighbor
      while (sum(duplicate >0) ) {

        # find 6 nearest host to each guest
        next_nn= lapply(1:num_neighbors, function(x) {
          nncross(upp_guest, host_unique, k = x)
        })

        ####
        dist_coords = lapply(1:num_neighbors, function(x) {
          coords(upp_guest) - coords(host_unique[next_nn[[x]][,2]])
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
      #dim(chosen_points)

      # get coordinates of all original points and their selected neighbors
      coords_multimer = rbind(coords(upp_guest)[,c(1,2, 3)],
                              data.frame(x = chosen_points$x, y = chosen_points$y, z = chosen_points$z))
      coords_multimer = subset(coords_multimer, z > min_thick &  z < max_thick )
      dim(coords_multimer)

      # make ppp and caluclate summary functions
      multimer_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
    }
    ### END CASE 1


    # case 2 UPP is 2d
    else if (spatstat.geom::is.ppp(upp)) {
      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2])
      ## label n/2 points as guest and the rest as host
      n_guest = n_guests
      n_total = npoints(upp_box)
      n_host = n_total - n_guest

      # determine the number of groups of molecules to be present (if multimers, then n_guest/2)
      n_groups = round(n_guest/ group_size,0)
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
      # extract guest points
      upp_guest = subset(upp_labeled, marks == "A")
      if (group_size == 1) {
        multimer_box = ppp(x = upp_guest$x, y = upp_guest$y, window = box_2d)
        return(multimer_box)
        next
      }
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")
      # find 6 nearest host to each guest
      nn_ind= lapply(1:num_neighbors, function(x) {
        nncross(upp_guest, upp_host, k = x)
      })
      # get nn distance x, y, and z values for each neighbor
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_guest) - coords(upp_host[nn_ind[[x]][,2]])
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


      ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two guests )
      # take the 2nd nearest neighbor
      while (sum(duplicate >0) ) {

        # find 6 nearest host to each guest
        next_nn= lapply(1:num_neighbors, function(x) {
          nncross(upp_guest, host_unique, k = x)
        })

        ####
        dist_coords = lapply(1:num_neighbors, function(x) {
          coords(upp_guest) - coords(host_unique[next_nn[[x]][,2]])
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
      coords_multimer = rbind(coords(upp_guest), data.frame(x = chosen_points$x, y = chosen_points$y))
      dim(coords_multimer)

      # make ppp and caluclate summary functions
      multimer_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
    }
    return(multimer_box)
    ## END CASE 2
  }

  # case 3: guest_pattern is 3d and upp is 3d
  else if (spatstat.geom::is.pp3(guest_pattern) || output == "pp3") {

    if (is.null(guest_pattern)) {
      box_3d = domain(upp)
    }

    else {
      box_3d = domain(guest_pattern)
      n_guests = npoints(guest_pattern)
    }

    box_area = spatstat.geom::volume(box_3d)

    if (is.na(max_thick)) {
      max_thick = box_3d$zrange[2]
    }
    if (is.na(min_thick)) {
      min_thick = box_3d$zrange[1]
    }



    # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
    upp_box = subset(upp_scaled, x >= box_3d$xrange[1] & x <= box_3d$xrange[2] &
                       y >= box_3d$yrange[1] & y <= box_3d$yrange[2] &
                       z >=min_thick - ztrim& z <= max_thick + ztrim)


    ## If using ztrim, then must adjust so that trimmed box has correct number of points
    # if not using ztrim, it will be zero and full/final will just be 1
    full_volume = abs(box_3d$xrange[2] - box_3d$xrange[1]) *
      abs(box_3d$yrange[2] - box_3d$yrange[1]) *
      (max_thick - min_thick + ztrim *2)
    final_volume = abs(box_3d$xrange[2] - box_3d$xrange[1]) *
      abs(box_3d$yrange[2] - box_3d$yrange[1]) *
      (max_thick - min_thick)


    # npoints will be scaled by volume ratios
    ## label n/2 points as guest and the rest as host
    n_guest = n_guests * full_volume/final_volume
    n_total = npoints(upp_box) #* full_volume/final_volume
    n_host = n_total - n_guest


    # determine the number of groups of molecules to be present (if multimers, then n_guest/2)
    n_groups = round(n_guest/ group_size,0)
    upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
    # extract guest points
    upp_guest = subset(upp_labeled, marks == "A")

    if (group_size == 1) {
      multimer_box = pp3(x = upp_guest$data$x, y = upp_guest$data$y,z = upp_guest$data$z, window = box_3d)
      ind = match(do.call("paste", coords(multimer_box)), do.call("paste", coords(upp_labeled)))
      marks(upp_labeled) = "H"
      marks(upp_labeled)[ind] = "G"
      marks(upp_labeled) = as.factor(marks(upp_labeled))
      return(upp_labeled)
      next
    }
    # extract host points
    upp_host = subset(upp_labeled, marks == "C")
    # find 6 nearest host to each guest
    nn_ind= lapply(1:num_neighbors, function(x) {
      nncross(upp_guest, upp_host, k = x)
    })
    # get nn distance x, y, and z values for each neighbor
    dist_coords = lapply(1:num_neighbors, function(x) {
      coords(upp_guest) - coords(upp_host[nn_ind[[x]][,2]])
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


    ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two guests )
    # take the 2nd nearest neighbor
    while (sum(duplicate >0) ) {

      # find 6 nearest host to each guest
      next_nn= lapply(1:num_neighbors, function(x) {
        nncross(upp_guest, host_unique, k = x)
      })

      ####
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_guest) - coords(host_unique[next_nn[[x]][,2]])
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
    coords_multimer = rbind(coords(upp_guest),
                            data.frame(x = chosen_points$x, y = chosen_points$y, z = chosen_points$z))
    coords_multimer = subset(coords_multimer, z > min_thick &  z < max_thick )
    #dim(coords_multimer)

    # make ppp and caluclate summary functions
    multimer_box = pp3(x = coords_multimer$x, y = coords_multimer$y, z = coords_multimer$z, window = box_3d)
    ind = match(do.call("paste", coords(multimer_box)), do.call("paste", coords(upp_box)))
    marks(upp_box) = "H"
    marks(upp_box)[ind] = "G"
    marks(upp_box) = as.factor(marks(upp_box))
    return(upp_box)
  }
}
