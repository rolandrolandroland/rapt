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
