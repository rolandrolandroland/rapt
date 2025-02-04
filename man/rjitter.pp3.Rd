% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rapt_extend.R
\name{rjitter.pp3}
\alias{rjitter.pp3}
\title{Random Perturbation of a Point Pattern}
\usage{
\method{rjitter}{pp3}(
  X,
  distribution = "uniform",
  xmin = 0,
  xmax = 0,
  ymin = 0,
  ymax = 0,
  zmin = 0,
  zmax = 0,
  xmean = 0,
  xsd = 0,
  ymean = 0,
  ysd = 0,
  zmean = 0,
  zsd = 0,
  retry = TRUE,
  giveup = 10000,
  ...,
  nsim = 1,
  drop = TRUE
)
}
\description{
Applies independent random displacements to each point in a point pattern.
Extends \code{\link[spatstat.geom]{rjitter}} to \code{\link[spatstat.geom]{pp3}}.
}
\seealso{
\code{\link[spatstat.geom]{rjitter}}, \code{\link{rjitter.ppp}}

Other spatstat extensions: 
\code{\link{G3cross}()},
\code{\link{G3est_nn}()},
\code{\link{G3multi}()},
\code{\link{K3scaled}()},
\code{\link{Tstat.pp3}()},
\code{\link{bdist.points}()},
\code{\link{marktable}()},
\code{\link{marktable.pp3}()},
\code{\link{quadratcount.pp3}()},
\code{\link{quadrats.pp3}()},
\code{\link{rPoissonCluster3}()},
\code{\link{rjitter}()},
\code{\link{rjitter.ppp}()},
\code{\link{rpoint3}()},
\code{\link{sample.pp3}()},
\code{\link{sample.ppp}()},
\code{\link{shift.pp3}()},
\code{\link{studpermu.test}()},
\code{\link{studpermu.test.pp3}()},
\code{\link{superimpose.pp3}()}
}
\concept{spatstat extensions}
