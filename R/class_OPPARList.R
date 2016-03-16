#' A S4 class for the output of OPPAR main function, opa.
#'
#' @slot profileMatrix A matrix of 0, -1 and 1 representing outlier genes in samples
#' @slot upper.quantile Numeric. The upper quantile cut-off for detection of outliers
#' @slot lower.quantile Numeric. The lower quantile cut-off for detection of outliers
#' @slot group A factor vector representing the group to which each sample belong
#' @slot .Data matrix
#' @export
setClass("OPPARList",
				 representation = representation(profileMatrix = "matrix", upper.quantile = "numeric",
				 																lower.quantile = "numeric", group = "factor"),
				 contains = "matrix")



#' @param object An object of type OPPARList
#' @describeIn OPPARList A show method for objects of class OPPARList
#' @export
setMethod("show","OPPARList", function(object){
	cat("Object of type OPPARList", sep="\n")
	cat(sprintf("Features: %d", dim(object@profileMatrix)[1]), sep = "\n")
	cat(sprintf("Samples: %d", dim(object@profileMatrix)[2]), sep = "\n")
	cat(sprintf("Upper quantile: %.2f", object@upper.quantile), sep = "\n" )
	cat(sprintf("Lower quantile: %.2f", object@lower.quantile), sep = "\n" )
	levels(object@group) <- c("0","1")
	cat("Groups:", sep = "\n")
	print(object@group)


})


#' @param x Object of type OPPARList
#' @param name Name of the slot to access
#' @describeIn OPPARList A method to extract slots in \code{OPPARList}
#' @export
setMethod("$", "OPPARList", function(x, name){
	eval(substitute(x@name))
	})



