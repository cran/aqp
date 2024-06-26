

.Gmedian.chunk <- function(i, id = '.chunk') {
  
  # iterate over chunks
  depth.list <- split(i, i[[id]])
  
  res <- lapply(depth.list, function(j) {
    
    # structure: group, chunk, top, bottom, var1, var2, var3, ...
    v <- 5:ncol(j)
    
    # NA not allowed, keep track and filter
    idx <- which(complete.cases(j[, v, drop = FALSE]))
    
    # catch conditions where there are no data
    # TODO: what is the min number of required records?
    if(length(idx) < 1)
      return(NULL)
    
    ## TODO: how do arguments to Gmedian() affect the results?
    # L1 median for non-NA records, variables of interest only
    G <- data.frame(Gmedian::Gmedian(X = j[idx, v, drop = FALSE]))
    
    # retain original names
    names(G) <- names(j)[v]
    
    # package group + top + bottom + L1 data
    d <- data.frame(group = j[1, 1], top = j[1, 3], bottom = j[1, 4], G, stringsAsFactors = FALSE)
    
    return(d)
  })
  
  # list -> DF
  res <- do.call('rbind', res)
  
  return(res)
}


## TODO:
# would it make more sense to select a single profile based on NCSP?
# how does variable scale affect Gmedian? 
# best practices for variable name compatibility (slab-style vs. original names)
# slab-style support for arbitrary depth intervals via basis argument
# most-likely horizon designation by chunk
# parallel Gmedian computation
# add principal components for viz
# document!

#' @title Create Representative Soil Profiles via L1 Estimator
#' 
#' @description The L1 estimator, or \href{https://en.wikipedia.org/wiki/Geometric_median}{geometric median}, is a multivariate generalization of the (univariate) \href{https://en.wikipedia.org/wiki/Median}{median} concept. This function performs a multivariate aggregation (via L1 estimator) according to a suite of ratio-scale soil properties. The L1 estimator is applied to soil profile data that have been sliced to a 1-depth-unit basis. Data should be well stratified by groups defined in `fm`, otherwise the L1 median may not make any sense.
#' 
#' See the \href{https://ncss-tech.github.io/AQP/aqp/L1-profiles.html}{L1 Profiles Tutorial} for additional examples.
#'
#' @note This function requires the `Gmedian` package.
#' 
#' @references Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. Bernoulli, 19, 18-43.
#'
#' @param x \code{SoilProfileCollection} object
#' 
#' @param fm formula, for example: `group ~ p1 + p2 + p3`, where "group" is a site-level grouping variable, and "p1", "p2", and "p3" are horizon level variables
#' 
#' @param basis positive integer, aggregation basis (e.g. 1 for 1-depth-unit intervals). Values other than 1 are not currently supported.
#' 
#' @param method soil depth evaluation method: "regex" for regular expression, "simple", or "constant". See details.
#' 
#' @param maxDepthRule maximum depth rule: "max" or "min" See details.
#' 
#' @param maxDepthConstant positive integer, maximum depth when \code{maxDepthRule = 'constant'}
#' 
#' @details See [this related tutorial](https://ncss-tech.github.io/AQP/aqp/L1-profiles.html) for additional examples. The `method`, `maxDepthRule`, and `maxDepthConstant` arguments set the maximum depth (over the entire collection) of analysis used to build "L1 profiles". The following rules are available:
#'   
#'   * `method = 'regex'` uses pattern matching on horizon designations (note that `hzdesgnname` metadata must be set with `hzdesgnname(x) <- 'columnname'`)
#'   
#'    * `method = 'simple'` uses `min` or `max` as applied to `x`, no accounting for non-soil horizons (e.g. Cr or R)
#'    
#'    * `method = 'constant'` uses a fixed depth value supplied by `maxDepthConstant`
#'
#' The `maxDepthRule` argument sets depth calculation constraint, applied to soil depths computed according to `method` (`min` or `max`).
#' 
#' @return a \code{SoilProfileCollection} object
#' @export
#'
L1_profiles <- function(x, fm, basis = 1, method = c('regex', 'simple', 'constant'), maxDepthRule = c('max', 'min'), maxDepthConstant = NULL) {
  
  # sanity check, need this for L1 median
  if(!requireNamespace('Gmedian')) {
    stop('package `Gmedian` is required', call. = FALSE)
  }
  
  # sanity checks: is this an SPC?
  if(! inherits(x, 'SoilProfileCollection')) {
    stop('`object` should be a SoilProfileCollection', call. = FALSE)
  }
  
  # extract components of the formula:
  g <- all.vars(update(fm, . ~ 0)) # left-hand side
  vars <- all.vars(update(fm, 0 ~ .)) # right-hand side
  
  # sanity check: do the variables specified in fm exist in the correct site / hz slots?
  if(
    ! all(vars %in% horizonNames(x)) |
    ! g %in% siteNames(x)
  ) {
    stop('`fm` must reference one or more horizon-level attributes and a single site-level attribute', call. = FALSE)
  }
  
  # arguments
  method <- match.arg(method)
  maxDepthRule <- match.arg(maxDepthRule)
  
  # multi-argument sanity
  if(method == 'constant' & !is.numeric(maxDepthConstant)) {
    stop('contant max depth must be specified by single numeric value', call. = FALSE)
  }
  
  # ensure that there is a horizon name defined when method = 'regex'
  if(method == 'regex'){
    if(hzdesgnname(x) == '') {
      stop('`x` must have a horizon designation defined to use `method = "regex"`, see `?hzdesgnname`', call. = FALSE)
    }
  }
  
  # SPC metadata
  hztb <- horizonDepths(x)
  
  # upper depth logic, usually 0
  collection.top <- min(x[[hztb[1]]], na.rm = TRUE)
  
  # lower depth logic
  collection.bottom <- switch(method,
         regex = {
           x.depths <- profileApply(x, estimateSoilDepth, name = hzdesgnname(x))
           switch(maxDepthRule, 
                  min = {
                    min(x.depths, na.rm = TRUE)
                  },
                  max = {
                    max(x.depths, na.rm = TRUE)
                  })
         },
         simple = {
           switch(maxDepthRule,
                  min = {
                    min(x)
                  },
                  max = {
                    max(x)
                  })
         },
         constant = {
           maxDepthConstant
         }
         
  )
   
  
  # create formula
  # this only needs the top/bottom depths
  slice.fm <- as.formula(sprintf("%s:%s ~ %s", collection.top, collection.bottom, paste(vars, collapse = ' + ')))
  
  # re-format for slice-wise L1 median
  s <- dice(x, fm = slice.fm, strict = TRUE, SPC = TRUE)
  
  # work on de-normalized data as a DF
  h <- as(s, 'data.frame')
  
  ## TODO: interpret basis as either slab thickness, or like slab.structure
  
  # determine which via length(basis)
  
  # assign slab IDs
  
  
  # simplest case: working with raw slices, use top depth
  if(basis == 1) {
    h$.chunk <- h[[hztb[1]]]
  } else {
    stop('sorry, depth basis other than `1` are not currently supported', call. = FALSE)
  }
  
  
  # retain only those columns we need
  h <- h[, c(g, '.chunk', hztb, vars)] 
  
  # iterate over groups
  group.list <- split(h, h[[g]])
  
  ## TODO: parallel
  agg <- lapply(group.list, .Gmedian.chunk)
  agg <- do.call('rbind', agg)
  
  # init SPC with new ID, top, bottom names
  depths(agg) <- group ~ top + bottom
  
  # transfer metadata
  agg <- .transfer.metadata.aqp(x, agg)
  
  # reset horizon designation for now
  hzdesgnname(agg) <- NULL
  
  return(agg)
}


