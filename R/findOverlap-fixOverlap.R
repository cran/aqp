#' Establish which elements within a vector of horizontal positions overlap beyond a given threshold
#'
#' @param x vector of relative horizontal positions, one for each profile
#' @return - `findOverlap` a vector of the same length as `x`, preserving rank-ordering and boundary conditions.
#' @note This is a very naive function and may fail to converge on a reasonable
#' solution. SANN would be a much more robust framework.
#' 
#' @rdname find-fix-overlap
#' 
#' @export
#'
findOverlap <- function(x, thresh) {
  # all pair-wise distance
  d <- dist(x)
  m <- as.matrix(d)
  # diagonal isn't used here
  diag(m) <- NA
  # find matrix elements
  idx <- which(m < thresh)
  # use upper-triangle indexes to find elements in original vector
  col.idx <- col(m)[idx]
  # done
  return(col.idx)
}

## 2019-07-16 | DEB
## fix overlap via random perterbation of affected elements
## this function is actually quite stupid as it can converge on bogus results
## scaled adjustments based on deviation from threshold distances would be better
## or, SANN style adjustments
##
## ideas:
## * cooling ala SANN -> smaller adjustments through time
## * perturbations should not increase overlap
## * debugging output
##
# x: vector of horizontal positions
# thresh: threshold at which overlap is a problem
# adj: adjustments are tested from runif(min=adj * -1, max=adj)
# min.x: left boundary condition
# max.x: right boundary condition
# maxIter: maximum number of iterations to attempt before collapsing to integer sequence
# trace: print diagnostics

#' Fix overlap via random perterbation of affected elements
#' @description This function is actually quite stupid as it can converge on bogus results
#' scaled adjustments based on deviation from threshold distances would be better
#  or, SANN style adjustments
#' @param x vector of horizontal positions
#' @param thresh horizontal threshold defining "overlap", must be < 1, ideal
#' values likely in \[0.3, 0.8]
#' @param adj adjustments are tested within \code{runif(min=adj * -1, max=adj)}
#' @param min.x left-side boundary condition
#' @param max.x right-side boundary condition
#' @param maxIter maximum number of iterations to attempt before giving up and
#' returning integer sequence
#' @param trace print diagnostics' 
#' @return - `findOverlap` a vector of the same length as `x`, preserving rank-ordering and boundary conditions.
#' @author D.E. Beaudette
#' @rdname find-fix-overlap
#' @export
#'
fixOverlap <- function(x, thresh=0.6, adj=0.2, min.x=0.8, max.x=length(x)+0.2, maxIter=1000, trace=FALSE) {
  
  # initial configuration
  ov <- findOverlap(x, thresh)
  
  # save original
  x.orig <- x
  
  # counter to prevent run-away while-loop
  i <- 1
  
  # short-circuit: only proceed if there is overlap
  if(length(ov) > 0) {
    
    # iterate...
    while(length(ov) > 0) {
      
      # fail-safe
      if(i > maxIter) {
        message('maximum number of iterations reached, using integer sequence')
        return(1:length(x))
      }
      
      # generate random perturbations to affected indices
      perturb <- runif(n = length(ov), min = adj * -1, max = adj)
      
      # attempt perturbation
      x.test <- x
      x.test[ov] <- x.test[ov] + perturb
      
      # enforce boundary conditions
      if(any(x.test < min.x) | any(x.test > max.x)) {
        # print('boundary condition')
        i <- i + 1
        next
      }
      
      # enforce rank ordering
      if(any(rank(x.orig) != rank(x.test))) {
        # print('rank violation')
        i <- i + 1
        next
      }
      
      ## this may be too strict: 85% -> 75% success rate if enabled
      # # perturbations should not increase number of affected positions
      # # check to be sure
      # len <- length(findOverlap(x.test, thresh))
      # # stats[[i]] <- len
      # if(len > length(ov)) {
      #   # print('wrong turn')
      #   i <- i + 1
      #   next
      # }
      
      ## alternative idea: perturbations should minimize overlap
      ## how to quantify?
      
      # apply perturbation to working copy
      x <- x.test
      
      # eval overlap and try again
      ov <- findOverlap(x, thresh)
      i <- i + 1
    }
  }
  
  if(trace) {
    message(sprintf("%s iterations", i))
  }
  
  
  return(x)
}

