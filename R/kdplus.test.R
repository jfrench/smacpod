#' Global test of clustering using difference in K functions
#' 
#' \code{kdplus.test} performs a global test of clustering for comparing cases and controls using the method of Diggle and Chetwynd (1991).  It relies on the difference in estimated K functions.
#' 
#' @param x A \code{kdenv} object from the \code{kdenv} function.  
#' 
#' @return A list providing the observed test statistic (\code{kdplus}) and the estimate p-value {\code{pvalue}}.
#' @author Joshua French
#' @import spatstat
#' @export
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.  Diggle, Peter J., and Amanda G. Chetwynd. "Second-order analysis of spatial clustering for inhomogeneous populations." Biometrics (1991): 1155-1163.
#' @examples 
#' data(grave)
#' env = kd.env(grave, nsim = 19)
#' kdplus.test(env)

kdplus.test = function(x)
{
  if(max(class(x) == "kdenv") < 1) stop("x must be an object from the kd.env function.")
  simfuns <- as.data.frame(attr(x, "simfuns"))
  sdkdhat = apply(simfuns, 1, sd) # estimated variance of kdest simulations
  # turn into matrix
  sdmat = matrix(sdkdhat, nrow = nrow(simfuns), ncol = ncol(simfuns))
  # estimate KD + for simulated data
  kdplussim = colSums(simfuns/sdmat, na.rm = TRUE)
  # determine proportion of simulated KD+ and observed KD+
  # greater than KD+
  kdplus = sum(x$obs/sdkdhat, na.rm = TRUE)
  p = (sum(kdplussim[-1] >= kdplus) + 1)/ncol(simfuns)
  print(paste("The p-value for the global test is", round(p, 3)))
  return(invisible(list(kdplus = kdplus, pvalue = p)))
}