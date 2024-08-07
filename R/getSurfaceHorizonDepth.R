#' Determine thickness of horizons (continuous from surface) matching a pattern
#'
#' @description Find the thickness of horizon designations, or any other character patterns, that are continuous from the soil surface (depth = 0 or shallowest depth in profile).
#'
#' @details The horizon designation to match is specified with the regular expression pattern 'pattern'. All horizons matching that pattern, that are continuous from the soil surface, count towards the depth / thickness value that is ultimately returned. For instance: horizon designations: A1-A2-A3-C-Ab , would return A3 bottom depth given \code{pattern = "^A[1-9]*$"}.
#'
#' `getSurfaceHorizonDepth` is used by `getPlowLayerDepth` for matching Ap horizons; and, it is used by `getMineralSoilSurfaceDepth` to find the thickness of O horizons in lieu of lab data.
#'
#' @param p a SoilProfileCollection
#' @param pattern a regular expression pattern to match for all horizons to be considered part of the "surface".
#' @param hzdesgn column name containing horizon designation. Default: `hzdesgnname(p, required = TRUE)`.
#' @param simplify logical. Return single profile results as vector (default: `TRUE`) or `data.frame` (`FALSE`)
#' @return a numeric value corresponding to the bottom depth of the last horizon matching 'pattern' that is contiguous with other matching horizons up to the soil surface. If `length(p) > 1` then a _data.frame_ containing profile ID, horizon ID, top or bottom depths, horizon designation and pattern.
#'
#' @author Andrew G. Brown
#'
#' @aliases getMineralSoilSurfaceDepth
#' @aliases getPlowLayerDepth
#'
#' @export 
#' @examples
#' library(aqp)
#' data(sp1, package = 'aqp')
#' depths(sp1) <- id ~ top + bottom
#' site(sp1) <- ~ group
#'
#' p <- sp1[1]
#' q <- sp1[2]
#'
#' # look at horizon designations in p and q
#' p$name
#' q$name
#'
#' # thickness of all surface horizons containing A
#' getSurfaceHorizonDepth(p, pattern = 'A', hzdesgn = 'name')
#'
#' # thickness of all surface horizons that start with A
#' getSurfaceHorizonDepth(p, pattern = '^A', hzdesgn = 'name')
#'
#' # thickness of all surface horizons that start with A, and the A is not followed by B
#' getSurfaceHorizonDepth(p, pattern = '^A[^B]*', hzdesgn = 'name')
#'
#' # thickness of all surface horizons that start with A
#' #  followed by a number from _2_ to 9 (returns ZERO)
#' getSurfaceHorizonDepth(p, pattern = '^A[2-9]*', hzdesgn = 'name')
#'
#' # getPlowLayerDepth matches first two horizons in fake Ap horizon data with "buried Ap"
#' p$aphorizons <- c("Ap1","Ap2","AB", rep('C', nrow(p) - 4), "Apb")
#' getPlowLayerDepth(p, hzdesgn = 'aphorizons')
#'
#' # getMineralSoilSurfaceDepthmatches first 3 horizons in fake O horizon data
#' p$ohorizons <- c("Oi1","Oi2","Oe", rep('C', nrow(p) - 4), "2C")
#' getMineralSoilSurfaceDepth(p, hzdesgn='ohorizons')
#'
#' # matches first Oi horizon with original horizon designations of pedon 2
#' getMineralSoilSurfaceDepth(q, hzdesgn='name')
#'
getSurfaceHorizonDepth <- function(p,
                                   pattern,
                                   hzdesgn = hzdesgnname(p, required = TRUE),
                                   simplify = TRUE) {
  
  
  if (is.null(hzdesgn) || !hzdesgn %in% horizonNames(p)) {
    stop("Horizon designation column (", hzdesgn, ") does not exist.")
  }
  
  hz <- data.table::as.data.table(horizons(p))
  depthz <- horizonDepths(p)
  
  .FIRST <- NULL
  .SD <- NULL
  
  shallowest.depth <- p[, , .FIRST][[depthz[1]]]
  
  .get_contiguous_surface_hz <- function(h, hzdesgn) {
    # identify horizons matching pattern
    match.idx <- grepl(h[[hzdesgn]], pattern = pattern)
    
    emptyres <- h[1, .SD, .SDcols = c(hzidname(p), depthz[2])] 
    emptyres[[idname(p)]] <- emptyres[[idname(p)]][1]
    emptyres[[hzidname(p)]] <- NA_character_
    emptyres[[depthz[2]]] <- as.double(0)
    emptyres$pattern <- pattern
    emptyres$hzdesgn <- hzdesgn
    
    # no match? return zero or shallowest top depth (the minimum depth)
    if (length(which(match.idx)) < 1) {
      return(emptyres)
    }
    
    # identify surface horizon
    mod.idx <- c(1, rep(0, length(match.idx) - 1))
    
    # identify where matches and surface horizon co-occur
    # matching horizons get a 1, matching surface horizon gets a 2,
    #   0s kick us out
    new.idx <- match.idx + mod.idx
    
    who.idx <- numeric(0)
    # we only have a matching surface if the first value is 2
    #  (meets both above crit)
    if (new.idx[1] == 2) {
      # convert that into logical to identify contiguous matches
      contig <- new.idx > 0 & new.idx <= 2
      
      # calculate difference between contiguous matches/nonmatches
      dcontig <- diff(as.integer(contig))
      
      # max depth is, at first, the bottom depth of last contiguous hz
      max.idx <- length(contig)
      
      if(length(max.idx)) {
        # if we have a negative change at any depth,
        #  we have a discontinuity
        
        # adjust max index (depth) accordingly
        if(any(dcontig < 0)) {
          # take first index of negative dcontig
          # add one to correct for indexing offset due to diff()
          max.idx <- which(dcontig < 0)[1] + 1
        }
        
        # return last value from contig
        # (last contiguous horizon with surface)
        who.idx <- rev(which(contig[1:max.idx]))[1]
      }
    }
    
    if (length(who.idx) == 0)
      return(emptyres)
    
    res <- h[who.idx, .SD, .SDcols = c(hzidname(p), depthz[2])] 
    res[[depthz[2]]] <- as.double(res[[depthz[2]]])
    res$pattern <- pattern
    res$hzdesgn <- hzdesgn
    res
  }
  
  hzids <- list(hz[[idname(p)]])
  names(hzids) <- idname(p)
  
  sitetemplate <- data.table::as.data.table(site(p))[,.SD, .SDcols = idname(p)]
  sitetemplate[[idname(p)]] <- as.character(sitetemplate[[idname(p)]])
  hzsub <- hz[, .get_contiguous_surface_hz(.SD, hzdesgn), by = hzids]
  hzsub[[idname(p)]] <- as.character(hzsub[[idname(p)]])
  out <- hzsub[sitetemplate, on = idname(p)]
  
  # not found
  notfound <- sum(shallowest.depth > out[[depthz[2]]], na.rm = TRUE)
  if (notfound > 0) {
    warning(sprintf("found %s profiles where pattern did not match and shallowest depth is greater than 0", notfound), call. = FALSE)
  }
  
  if(length(p) == 1 && simplify) {
    return(out[[depthz[2]]])
  }
  as.data.frame(out)
}

#' @rdname getSurfaceHorizonDepth
#' @export
getMineralSoilSurfaceDepth <-  function(p, hzdesgn = hzdesgnname(p, required = TRUE), pattern = "O", simplify = TRUE) {
  
  if (is.null(hzdesgn) || !hzdesgn %in% horizonNames(p)) {
    stop("Horizon designation column (", hzdesgn, ") does not exist.")
  }

  # assumes OSM is given O horizon designation;
  # TODO: add support for lab-sampled organic measurements
  #       keep organic horizons with andic soil properties
  return(getSurfaceHorizonDepth(p, hzdesgn = hzdesgn, pattern = pattern, simplify = simplify))
}

#' @rdname getSurfaceHorizonDepth
#' @export
getPlowLayerDepth <- function(p, hzdesgn = hzdesgnname(p, required = TRUE), pattern = "^Ap[^b]*", simplify = TRUE) {
  return(getSurfaceHorizonDepth(p, hzdesgn = hzdesgn, pattern = pattern, simplify = simplify))
}
