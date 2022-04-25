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
