
#' @title Probabilistic Estimation of Soil Depth within Groups
#'
#' @description Estimate the most-likely depth to contact within a collection of soil profiles. Consider `getSoilDepthClass` followed by group-wise percentile estimation as a faster alternative.
#'
#' @param x a \code{SoilProfileCollection} object
#' @param groups the name of a site-level attribute that defines groups of profiles within a collection
#' @param crit.prob probability cutoff used to determine where the most likely depth to contact will be, e.g. 0.9 translates to 90% of profiles are shallower than this depth
#' @param name horizon-level attribute where horizon designation is stored, defaults to `hzdesgnname(x)`
#' @param p a REGEX pattern that matches non-soil genetic horizons
#' @param ... additional arguments to \code{slab}
#'
#' @details This function computes a probability-based estimate of soil depth by group. If no grouping variable exists, a dummy value can be used to compute a single estimate. The \code{crit.prob} argument sets the critical probability (e.g. 0.9) at which soil depth within a group of profiles is determined. For example, a \code{crit.prob} of 0.95 might result in an estimated soil depth (e.g. 120cm) where 95% of the profiles (by group) had depths that were less than or equal to 120cm.
#'
#' @return A `data.frame` is returned, with as many rows as there are unique group labels, as specified in `groups`.
#'
#' @author D.E. Beaudette
#'
#' @seealso [estimateSoilDepth()] [slab()]
#'
#' @export
#'
#' @examples
#'
#' data(sp1)
#' depths(sp1) <- id ~ top + bottom
#' site(sp1) <- ~ group
#'
#' # set horizon designation in SPC
#' hzdesgnname(sp1) <- 'name'
#'
#' aggregateSoilDepth(sp1, 'group', crit.prob = 0.9)
#'
aggregateSoilDepth <- function(x, groups, crit.prob = 0.9, name = hzdesgnname(x), p = 'Cr|R|Cd', ...) {

  # sanity checks:
  # * does the group label variable exist in @site?
  if(! groups %in% siteNames(x))
    stop('`groups` must specify a site-level attribute')

  # * does the horizon name variable exist in @horizons?
  if(! name %in% horizonNames(x))
    stop('`name` must specify a horizon-level attribute, containing the horizon designation')

  # mark soil vs. non-soil horizons
  x$soil.flag <- rep('soil', times = nrow(x))
  x$soil.flag[grep(p, horizons(x)[[name]])] <- 'not-soil'

  # convert to factor
  x$soil.flag <- factor(x$soil.flag, levels = c('soil', 'not-soil'))

  ## TODO: this is the slowest step
  # compute slice-wise probabilities that sum to contributing fraction
  fm <- as.formula(paste0(groups, ' ~ soil.flag'))
  a.ml.soil <- slab(x, fm, cpm = 2, ...)

  ## note: this is important when aggregating over several depth classes
  # compute adjusted soil flag probability: cf * Pr(soil)
  a.ml.soil$adj.soil.flag <- with(a.ml.soil, contributing_fraction * soil)

  # convert probability of contact into probability of soil
  crit.prob <- 1 - crit.prob

  # extract depth at specific probability of contact
  ss <- split(a.ml.soil, a.ml.soil[[groups]])
  depth.prob <- lapply(ss, function(i) {

    soil.top <- 0
    soil.bottom <- max(i$top[which(i$adj.soil.flag >= crit.prob)])

    res <- data.frame(
      .id = i[[groups]][1],
      soil.top,
      soil.bottom,
      stringsAsFactors = FALSE
    )

    names(res)[1] <- groups
    return(res)
  })

  # flatten list -> data.frame, reset rownames
  depth.prob <- do.call('rbind', depth.prob)
  row.names(depth.prob) <- NULL

  return(depth.prob)
}

