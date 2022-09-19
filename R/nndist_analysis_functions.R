#' nndist plotter and saver
#' @param image Point pattern for analysis (ppp or pp3)
#' @param image2 If a point pattern is given for image2, then \code{\link[spatstat]{nncross}} will be used
#' @param nn_vec vector of nearest neighbors to calculate distance to
#' @param location location to save image(s) if `output` = `save`
#' @param image_name Name to save image as.  The nearest neighbor number + `nn.png` will be added to it
#' @param xlim limits of x axis of histogram
#' @param ylim limits of y axis of histogram
#' @param binSize histogram bin size
#' @param freq argument for  \code{\link{hist}} function. logical;
#' if TRUE, the histogram graphic is a representation of frequencies, the counts component of the result;
#' if FALSE, probability densities, component density, are plotted (so that the histogram has a total area of one).
#'  Defaults to TRUE if and only if breaks are equidistant (and probability is not specified).
#' @param image_width Width of image to be saved
#' @param image_height Height of image to be saved
#'
#' @description
#' Plot and Save Histograms of nearest neighbor distances (\code{\link[spatstat]{nndist}})
#'
#' @details
#' Takes a point pattern as input (or two if using multitype \code{\link[spatstat]{nncross}}) and returns the histograms of all nearest neighbors
#' in the nn_vec argument.  If `output` is set to "save", then it will also save the image.
#'
#'

#'
#' @return Describe output
#' @export
nndist_func_plotter = function(image, image2 = NULL, nn_vec = 1:5, location = getwd(),
                               image_name = "image", output = "save",
                               xlim = c(0, 10), ylim = c(0, 0.5),
                               binSize = 0.5, freq = FALSE,
                               image_width = 600, image_height = 500,
                               ...) {
  # if image2 is null, then use nndist
  if (is.null(image2)) {
    distances = nndist(image, k = nn_vec)
  }
  else {
    # if there is an input pattern for image 2, then use nncross(image, image2)
    distances = nncross(image, image2, k = nn_vec)[,1:length(nn_vec)]
  }
  # define empty list to store output
  ls <- vector(mode = "list", length = length(nn_vec))
  # if you want to save, rather than just plot
  if (output == "save") {

    for(nn in 1:length(nn_vec)) {
      # the histogram will have range from the smallest distance, rounded down, to the largest bin, rounded up
      # note: this inherently requires the distances to be at least on the order of 10^0
      start = floor(min(distances[,nn]))
      stop = ceiling(max(distances[,nn]))
      # save at location `location` with name `image_name`_`nn`_nn.png
      png(file= paste(location, "/",image_name, "_",nn, "nn.png", sep = ""),
          width=image_width, height=image_height)
      hist(distances[,nn], xlim =xlim, ylim = ylim, breaks = seq(start, stop, binSize),
           freq = freq, main = paste(image_name, "Distance for NN", nn_vec[nn]), xlab = "Distance")
      dev.off()
      # store in output for plotting
      ls[[nn]] =hist(distances[,nn], xlim = xlim, ylim = ylim, breaks = seq(start, stop, binSize),
                     freq = freq, main = paste(image_name, "Distance for NN", nn_vec[nn]), xlab = "Distance")
    }
  }
  # if only plotting, do same as above but skip saving step
  else if (output == "plot") {
    for(nn in 1:length(nn_vec)) {

      start = floor(min(distances[,nn]))
      stop = ceiling(max(distances[,nn]))
      ls[[nn]] =hist(distances[,nn], xlim = xlim, ylim = ylim, breaks = seq(start, stop, binSize),
                     freq = freq, main = paste(image_name, "Distance for NN", nn_vec[nn]), xlab = "Distance")
    }
  }
  else {
    print("error: output should either be 'plot' or 'save'")
  }
  return(ls)
}

#' Plot and Save Histograms of nndist function
#'
#'  @param image Image of type ppp or pp3 to calculate distances from
#'  @param image2 will be NULL (default) if using \code{\link[spatstat]{nndist}} function on `image`.  Otherwise, set as same class as `image` (ppp or pp3) to
#'  calculate distances to using (\code{\link[spatstat]{nncross}})
#'  @param value determines the output.  If = "percent", then rows are NN number, columns are distance, and value is percent.  If = "distance",
#'  rows are NN number, columns are percent, and `value` is distance
#'  @param rads_vec vector of distances to calculate percents at. Only used if `value` = "percent"
#'  @param nn_vec vector of integers values to determine which nearest neighbors to calculate
#'  @param perc_vec vector of percent values between 0 and 1 (inclusive) to calculate distances for.  Only used if `value` = "distance"
#'  @param output set to "plot" (default) if only plotting.  Set to "save" if plotting and saving
#'  @param location location to save image to if `output = save`.  default is working directory.
#'  @param image_name Name to save image as.  Will have `_NN_dist_perc.png` added to it if `value` = "percent" or `_NN_perc_dist.png` if `value` = "distance"
#'  @param unit distance unit name to be added to plot name
#'  @param round_to number of decimal places to round values to
#'
#' @description
#' calculate percent of each nearest neighbor (NN) in `nn_vec` that are less than x `unit` in `rads_vec` apart if `value` = "percent" or calculate the distance
#' at which each percentage in `perc_vec` has each NN in `nn_vec` if `value` = "distance."
#'
#'
#' @details
#'  add details
#'

#'
#' @return returns at matrix. If `value` = "percent", then rows are NN number, columns are distance, and value is percent.  If `value` = "distance",
#'  rows are NN number, columns are percent, and `value` is distance.  If `output` = "save" then it will also save the image at location `location`
#' @export
nn_dist_perc_table = function(pattern, pattern2 = NULL,
                              window = NULL, drop_isolated = FALSE,
                              value = "percent",
                              rads_vec = 1:25, nn_vec = 1:5,
                              perc_vec = seq(0, 1, 0.1), output = "plot",
                              location = getwd(),
                              pattern_name = "pattern",
                              unit = "nm", round_to = 2, ...) {
  distances = nndist_subset(pattern, pattern2, window, drop_isolated, k = nn_vec, output = "matrix")
  distances = distances[complete.cases(distances),]
  # if the desired matrix data type is percent
  if (value == "percent") {
    # calculate the percent of values that are smaller than each rads_vec value for each nn_vec value
    out = sapply(rads_vec, function(rad) {
      sapply(1:length(nn_vec), function(nn) {
        sum(distances[,nn] < rad)/length(distances[,nn])
      })
    })
    colnames(out) = sapply(rads_vec, function(x) { paste(x, unit)})
    rownames(out) =  sapply(nn_vec, function(x) { paste(x, "NN")})
    out = round(out, round_to)

    # if only plotting and not saving, then just return it
    if (output == "plot") {
      return(out)
    }
    # if saving, then make prettier using kable then save
    else if (output == "save") {
      print("Save")
      out %>% kable(caption = "Percent with all NN within Distance") %>%
        kable_styling(latex_options = c("striped")) %>%
        save_kable(paste(location, "/", pattern_name, "_NN_dist_perc.png", sep = ""))
      return(out)
    }
    else {
      print("Error: output must be either 'print' or 'save'")
    }
  }
  # if matrix values are distances
  else if (value == "distance") {
    # calulate the distances at which each percentage in perc_vec has each NN in nn_vec
    out = apply(distances[,nn_vec], 2, quantile, probs = perc_vec)
    colnames(out) =  sapply(nn_vec, function(x) { paste(x, "NN")})
    out = round(out,round_to)
    out = t(out)

    # plot or save and plot, same as above
    if (output == "plot") {
      return(out)
    }
    else if (output == "save") {
      print("Save")
      out %>% kable(caption = "Distance at Which X% of points have Y NN's") %>%
        kable_styling(latex_options = c("striped")) %>%
        save_kable(paste(location, "/", pattern_name, "_NN_perc_dist.png", sep = ""))
      return(out)
    }
    else {
      print("Error: output must be either 'print' or 'save'")
    }
  }
  else {
    print("`value` must be either `distance` or `percent`")
  }
}

#' Fraction of all neighbors within a distance that are same type or fraction of all nth nearest neighbors that are same type
#' @param pattern Image of type ppp or pp3
#' @param i mark for point type to be used
#' @param value Determines output.  If = "NN" then returns fraction of all nearest neighbors in nn_vec that are same type.  If = "distance" then returns fraction of
#'  all points within distances in rads_vec that are same type
#' @param rads_vec vector of distances to calculate
#' @param nn_vec vector of integers values of nearest neighbors
#' @param output set to "plot" (default) if only plotting.  Set to "save" if plotting and saving
#' @param location location to save image to if `output = save`.  default is working directory.
#' @param image_name Name to save image as.  Will have `_NN_dist_perc.png` added to it if `value` = "percent" or `_NN_perc_dist.png` if `value` = "distance"
#' @param unit distance unit name to be added to output name
#' @param round_to number of decimal places to round values to
#'
#' @description Calculate either the average fraction of points in object `image` at distances `rads_vec` from each point of type `i` that are type `i` (`value` = "distance") or the average fraction
#' of each nearest neighbor in `nn_vec` of each point of type `i` that is type `i` `value` = "NN").
#'
#' @return vector of fractions
#' @export
percent_neighbors = function(pattern, i,
                             window = NULL, drop_isolated = FALSE,
                             value = "NN",
                             rads_vec = 1:25, nn_vec = 1:5,
                            output = "plot",
                             location = getwd(),
                             image_name = "image",
                             unit = "nm", round_to = 3, ...) {
  # determine which output to use
  if (value == "NN") {

        # find distance from each point of type `i` to NN
    dist = nndist_subset(pattern[pattern$data$marks ==i], pattern, k = nn_vec, window = window,
                         drop_isolated = drop_isolated, output = "matrix")
    dist = dist[complete.cases(dist),]
    # fraction of NN's that is type `i`
    percent_neighbors = apply(dist[,(length(nn_vec)+1):(length(nn_vec)*2)], 2, function(x) {
      type =pattern$data$marks[x]
      total =sum(type ==i)
      total/ nrow(dist)
    })

    # round and name output
    out =round(percent_neighbors, round_to)
    names(out) = paste(nn_vec, "NN")
    name = "Fraction of each NN that is same type"
  }
  # determine if using `distance` output
  else if(value == "distance")
  {
    # all pairs of type `i` that are within each distance in `rads_vec`
    pairs_i = lapply(rads_vec, function(x) {
      closepairs(pattern[pattern$data$marks ==i], x)
    })
    # total number of pairs per distance
    num_pairs_i = sapply(pairs_i, function(x) {length(x[[1]])})

    # pairs of points that involves at least one point of type `i` within each distance in `rads_vec`
    pairs_i_all = lapply(rads_vec, function(x) {
      crosspairs(pattern[pattern$data$marks ==i], pattern, x)
    })

    # total number of pairs_i_all per distance
    num_pairs_i_all = sapply(pairs_i_all, function(x) {length(x[[1]])})
    # must subtract off all of the pairs that are self - self
    num_pairs_i_all = num_pairs_i_all - sum(pattern$data$marks ==i)
    # fraction that are type i to type -i
    out = num_pairs_i /num_pairs_i_all
    names(out) = paste(rads_vec, unit)
    out = round(out, round_to)
    name = "Fraction of all NN's within Distance that are same type"
    #return(out)
  }
  else {
    print("Value must be either 'distance' or 'NN' ")
  }

  if (output == "plot") {
    return(out)
  }

  # save the plot
  else if (output == "save") {
    print("Save")
    out %>% t()%>% kable(caption = name) %>%
      kable_styling(latex_options = c("striped")) %>%
      save_kable(paste(location, "/", image_name, "_NN_type_perc.png", sep = ""))
    return(out)
  }
  else {
    print("output must be either 'plot' or 'save'")
  }
}

#' Average data frames
#' @description
#' add description
#'
#' @details
#'  add details
#'
#'  @param image
#'  @param image2
#'
#'
#' @return Describe output
#' @export
data_frame_list_averager = function(data, output = "plot",
                                    location = getwd(),
                                    image_name = "image") {
  mean_df = data.frame(matrix(data = NA, nrow = nrow(data[[1]]), ncol =ncol(data[[1]])))
  sd_df  = data.frame(matrix(data = NA, nrow = nrow(data[[1]]), ncol = ncol(data[[1]])))
  for (col in 1:ncol(mean_df)) {
    for (row in 1:nrow(mean_df)) {
      # get all values into one vector
      vals = sapply(data, function(x) {
        x[row, col]
      })
      mean_df[row, col] = round(mean(vals), 2)
      sd_df[row, col] = round(sd(vals), 2)
    }
    out = list("mean_df" = mean_df, "sd_df" = sd_df)
    colnames(out[[1]]) = colnames(data[[1]])
    rownames(out[[1]]) =  rownames(data[[2]])
    colnames(out[[2]]) = colnames(data[[1]])
    rownames(out[[2]]) =  rownames(data[[2]])
  }
  if (output == "plot") {
    return(out)
  }
  else if (output == "save") {
    print("yaa")
    out[[1]] %>% kable(caption = "Mean Percent with all NN within Distance") %>%
      kable_styling(latex_options = c("striped")) %>%
      save_kable(paste(location, "/", image_name, "_mean_perc_table.png", sep = ""))
    out[[2]] %>% kable(caption = "Standard Deviation of Percent with all NN within Distance") %>%
      kable_styling(latex_options = c("striped")) %>%
      save_kable(paste(location, "/", image_name, "_sd_perc_table.png", sep = ""))
    return(out)
  }
  else {
    print("Error: output must be either 'print' or 'save'")
  }
}
