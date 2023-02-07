# I removed this so that rjitter.pp3 actually workds like rjitter.ppp
rjitter.pp3 <- function(X, domain = box3()) {
  spatstat.geom::verifyclass(X, "pp3")
  nX <- spatstat.geom::npoints(X)
  if (nX == 0) {
    return(X)
  }
  W <- X$domain
  D <- spatstat.random::runifpoint3(nX, domain = domain)
  xnew <- X$data$x + D$data$x
  ynew <- X$data$y + D$data$y
  znew <- X$data$z + D$data$z
  new <- pp3(xnew, ynew, znew, W)
  ok <- subset(new,
               subset =
                 (x > W$xrange[1] & x < W$xrange[2]) &
                 (y > W$yrange[1] & y < W$yrange[2]) &
                 (z > W$zrange[1] & z < W$zrange[2])
  )
  return(ok)
}

rjitter.pp3 <- function(X, distribution = "uniform",
                        xmin = 0, xmax = 0,
                        ymin = 0, ymax = 0,
                        zmin = 0, zmax = 0,
                        xmean = 0, xsd = 0,
                        ymean = 0, ysd = 0,
                        zmean = 0, zsd = 0,
                        retry = TRUE, giveup = 10000, ...,
                        nsim = 1, drop = TRUE
                        ) {
  verifyclass(X, "pp3")
  if (!missing(nsim)) {
    check.1.integer(nsim)
    stopifnot(nsim >= 0)
  }
  nX <- npoints(X)
  W <- domain(X)
  if (nX == 0) {
    result <- rep(list(X), nsim)
    result <- simulationresult(result, nsim, drop)
    return(result)
  }

  #'
  result <- vector(mode = "list", length = nsim)
  for (isim in seq_len(nsim)) {
    if (!retry) {

      ## points outside window are lost
      # uniform distribution between 0 and the radius
      if (distribution == "uniform") {
        dX = runif(nX, xmin, xmax)
        dY = runif(nX, ymin, ymax)
        dZ = runif(nX, zmin, zmax)

      }
      else if (distribution == "normal") {
        dX = rnorm(nX, xmean, xsd)
        dY = rnorm(nX, ymean, ysd)
        dZ = rnorm(nX, zmean, zsd)
      }
      else {
        stop("Sorry, currently only `uniform` and `normal` distrubutions are supported")
      }

      xnew <- X$data$x + dX
      ynew <- X$data$y + dY
      znew <- X$data$z + dZ
      ok <- inside.boxx(xnew, ynew, znew, w = W)
      result[[isim]] <- pp3(xnew[ok], ynew[ok], znew[ok], W)
    } else {
      ## retry = TRUE: condition on points being inside window
      undone <- rep.int(TRUE, nX)
      triesleft <- giveup
      Xshift <- X
      while (any(undone)) {
        triesleft <- triesleft - 1
        if (triesleft <= 0) {
          break
        }
        Y <- Xshift[undone]
        nY <- npoints(Y)
        print(nY)
        #RY <- if (sameradius) radius else radius[undone]

        ## points outside window are lost
        # uniform distribution between 0 and the radius
        if (distribution == "uniform") {
          dX = runif(nY, xmin, xmax)
          dY = runif(nY, ymin, ymax)
          dZ = runif(nY, zmin, zmax)

        }
        else if (distribution == "normal") {
          dX = rnorm(nY, xmean, xsd)
          dY = rnorm(nY, ymean, ysd)
          dZ = rnorm(nY, zmean, zsd)
        }
        else {
          stop("Sorry, currently only `uniform` and `normal` distrubutions are supported")
        }

        ## if there is only one point, then data structure is screwed up
        if (nY == 1) {
          xnew <- Y$data$x[[1]] + dX
          ynew <- Y$data$y[[1]] + dY
          znew <- Y$data$z[[1]] + dZ
          ok <- inside.boxx(xnew, ynew, znew, w = W)
        } else {
          xnew <- Y$data$x + dX
          ynew <- Y$data$y + dY
          znew <- Y$data$z + dZ
          ok <- inside.boxx(xnew, ynew, znew, w = W)
        }

        # result[[isim]] <- pp3(xnew[ok], ynew[ok], znew[ok], W)



        if (any(ok)) {
          changed <- which(undone)[ok]
          Xshift$data$x[changed] <- xnew[ok]
          Xshift$data$y[changed] <- ynew[ok]
          Xshift$data$z[changed] <- znew[ok]
          undone[changed] <- FALSE
        }
      }
      result[[isim]] <- Xshift
    }
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}


rjitter_dev <- function(X, radius, retry = TRUE, giveup = 10000, ...,
                        nsim = 1, drop = TRUE, dims = c("x", "y", "z"),
                        distribution = "uniform", mean = 0, sd = 1) {
  verifyclass(X, "pp3")
  if (!missing(nsim)) {
    check.1.integer(nsim)
    stopifnot(nsim >= 0)
  }
  nX <- npoints(X)
  W <- domain(X)
  if (nX == 0) {
    result <- rep(list(X), nsim)
    result <- simulationresult(result, nsim, drop)
    return(result)
  }
  if (missing(radius) || is.null(radius)) {
    ## Stoyan rule of thumb
    bws <- 0.15 / sqrt(5 * nX / volume(W))
    radius <- min(bws, shortside(W)) # had to get rid of "Frame" function call
    sameradius <- TRUE
  } else {
    ## either one radius, or a vector of radii
    check.nvector(radius, nX, oneok = TRUE, vname = "radius")
    sameradius <- (length(radius) == 1)
  }
  #'
  result <- vector(mode = "list", length = nsim)
  for (isim in seq_len(nsim)) {
    if (!retry) {

      ## points outside window are lost
      # uniform distribution between 0 and the radius
      if (distribution == "uniform") {
        rD <- radius * sqrt(runif(nX))
      }
      else if (distribution == "normal") {
        rD <- radius * sqrt(abs(rnorm(nX, mean, sd)))
      }
      else {
        stop("Sorry, currently only `uniform` and `normal` distrubutions are supported")
      }
      if ("z" %in% dims) {
        thetaD <- runif(nX, max = pi)
      }
      else {
        thetaD <- pi/2
      }
      if ("x" %in% dims) {
        if ("y" %in% dims) {
          phiD <- runif(nX, max = 2 * pi)

        }
        else {
          phiD <- 0
        }
      }
      # will only run if y but not x
      else if ("y" %in% dims) {
        phiD = pi/2}

      ## points outside window are lost
      #rD <- radius * sqrt(runif(nX))


      #thetaD <- runif(nX, max = pi)
      #phiD <- runif(nX, max = 2 * pi)
      xnew <- X$data$x + rD * sin(thetaD) * cos(phiD)
      ynew <- X$data$y + rD * sin(thetaD) * sin(phiD)
      znew <- X$data$z + rD * cos(thetaD)
      ok <- inside.boxx(xnew, ynew, znew, w = W)
      result[[isim]] <- pp3(xnew[ok], ynew[ok], znew[ok], W)
    } else {
      ## retry = TRUE: condition on points being inside window
      undone <- rep.int(TRUE, nX)
      triesleft <- giveup
      Xshift <- X
      while (any(undone)) {
        triesleft <- triesleft - 1
        if (triesleft <= 0) {
          break
        }
        Y <- Xshift[undone]
        nY <- npoints(Y)
        RY <- if (sameradius) radius else radius[undone]

        if (distribution == "uniform") {
          rD <- radius * sqrt(runif(nY))
        }
        else if (distribution == "normal") {
          rD <- radius * sqrt(abs(rnorm(nY, mean, sd)))
        }
        else {
          stop("Sorry, currently only `uniform` and `normal` distrubutions are supported")
        }
        if ("z" %in% dims) {
          thetaD <- runif(nY, max = pi)
        }
        else {
          thetaD <- pi/2
        }
        if ("x" %in% dims) {
          if ("y" %in% dims) {
            phiD <- runif(nY, max = 2 * pi)

          }
          else {
            phiD <- 0
          }
        }
        # will only run if y but not x
        else if ("y" %in% dims) {
          phiD = pi/2}
        # rD <- RY * sqrt(runif(nY))
        #thetaD <- runif(nY, max = pi)
        #phiD <- runif(nY, max = 2 * pi)

        ## if there is only one point, then data structure is screwed up
        if (nY == 1) {
          xnew <- Y$data$x[[1]] + rD * sin(thetaD) * cos(phiD)
          ynew <- Y$data$y[[1]] + rD * sin(thetaD) * sin(phiD)
          znew <- Y$data$z[[1]] + rD * cos(thetaD)
          ok <- inside.boxx(xnew, ynew, znew, w = W)
        } else {
          xnew <- Y$data$x + rD * sin(thetaD) * cos(phiD)
          ynew <- Y$data$y + rD * sin(thetaD) * sin(phiD)
          znew <- Y$data$z + rD * cos(thetaD)
          ok <- inside.boxx(xnew, ynew, znew, w = W)
        }

        # result[[isim]] <- pp3(xnew[ok], ynew[ok], znew[ok], W)



        if (any(ok)) {
          changed <- which(undone)[ok]
          Xshift$data$x[changed] <- xnew[ok]
          Xshift$data$y[changed] <- ynew[ok]
          Xshift$data$z[changed] <- znew[ok]
          undone[changed] <- FALSE
        }
      }
      result[[isim]] <- Xshift
    }
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}


