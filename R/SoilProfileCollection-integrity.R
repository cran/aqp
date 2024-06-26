#' Quickly assess relative state of site and horizon slots
#'
#' @description
#' 
#' Determine "state" of SoilProfileCollection before or after major modifications of site or horizon slot contents. 
#' 
#' Two logical checks are performed on the site and horizon tables, and a third element \code{valid} returns \code{TRUE} when both checks are \code{TRUE}.
#'
#' Check 1: Same number of sites in site as number of sites in horizons. No intermingling of IDs, no orphan horizons, no sites without horizons (for now)
#' 
#' Check 2: Site IDs match coalesced profile ID from horizons. Ensures the same *relative* ordering, but horizons still may be out of order within profiles
#' 
#' @param object A SoilProfileCollection 
#' @return data.frame
#' @docType methods
#' @author Andrew G. Brown
#' @export
#' @examples 
#' data(sp5)
#' 
#' spc_in_sync(sp5)
#' 
spc_in_sync <- function(object) {
  
  # get site and horizon slot contents
  s <- object@site
  h <- object@horizons
  
  # object profile ID name
  oid <- idname(object)
  hzidnm <- hzidname(object)
  
  # profile ID from site
  sid <- as.character(s[[oid]])
  
  # profile ID from horizon
  hid <- as.character(h[[oid]])
  hzid <- as.character(h[[hzidnm]])

  if (length(hid) == 0 & length(hzid)  == 0)  {
    # this is an empty soil profile collection
    return(data.frame(siteDepth = TRUE,
                    relativeOrder = TRUE,
                    depthOrder = TRUE,
                    valid = TRUE))
  }
  
  # coalesced horizon IDs 
  # identifies intermingling of profiles within horizon
  cohid <- rle(hid)$values
  
  # if cohid is longer than sid, horizons from different profiles
  # are mixed or IDs have been corrupted (by e.g. direct edit)
  # in less-common case, a site has been removed or duplicated
  
  # check 1: site order is reflected in horizons 
  #          no intermingling of sites within horizon table
  #          though, they _may_ be out of order within profiles
  one <- (length(sid) == length(cohid))
  
  # check 2: site IDs match coalesced profile ID from horizon
  #          this ensures the same _relative_ ordering, but
  #          horizons still may be out of order within profiles
  #          use checkHzDepthLogic -- as only one problematic pedon will break SPC depth order
  two <- ifelse(one, all(sid == cohid), FALSE)
  
  return(data.frame(nSites = one,
                    relativeOrder = two,
                    valid = all(one, two)))
}

# if (!isGeneric('reorderHorizons'))
  setGeneric('reorderHorizons', 
             function(object, target.order = NULL)
    standardGeneric('reorderHorizons'))

#' Re-order corrupted horizon data
#' 
#' @name reorderHorizons
#' @description This is a method made available primarily to repair horizon data that have been corrupted relative to their order at time of SoilProfileCollection construction. 
#' 
#' There is an option to specify the target order, but this will not update the corresponding metadata entry tracking the original order. Use this functionality at your own risk.
#' 
#' @param object A SoilProfileCollection
#' @param target.order A numeric vector of equal length to \code{object}. Default value is \code{NULL} which restores the internal order of the collection. 
#' @aliases reorderHorizons,SoilProfileCollection-method
#' @return SoilProfileCollection
#' @docType methods
#' @rdname reorderHorizons
#' @export
setMethod('reorderHorizons',
          signature(object = 'SoilProfileCollection', target.order = "ANY"),
          function(object, target.order = NULL) {
            
            h <- object@horizons
            
            stopifnot(!is.null(target.order))
            
            current.order <- match(target.order,
                                   order(as.character(h[[idname(object)]]),
                                         h[[horizonDepths(object)[1]]]))
            
            h <- .as.data.frame.aqp(h[current.order,], 
                                          aqp_df_class(object))
            object@horizons <- h
            return(object)
          })
