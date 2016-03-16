#' @include class_OPPARList.R class_nuchar.R
NULL

#' Retrieving outlier genes in samples
#'
#' Returns a list of outlier genes for each given sample
#'
#' @param profileMatrix  A matrix of 0, -1 and 1 representing outlier genes in samples.
#'   Also an object of type \code{OPPARList}.
#' @param sample.name A character vector containing one or more sample names, or a numeric vector
#'   containing sample indices.
#' @return A list of lists. For each given sample, the fuction return
#'    up-regulated and down-regulated outlier genes.
#' @export
#' @examples
#' data(GSE46141g)
#' library(Biobase)
#' eset <- g[[1]]
#' group <- sapply(pData(eset)$source_name_ch1, function(x){ ifelse(x == "breast",0,1)})
#' group <- factor(group)
#' bcm.opa <- opa(eset,group=group)
#' # Outlier profile for sample "GSM1124929"
#' getSampleOutlier(bcm.opa, "GSM1124929")
#'
#' # Also can use sample index, instead of sample name
#' getSampleOutlier(bcm.opa, 11)
#'
#' # A vector of sample names/indices are accepted as well
#' getSampleOutlier(bcm.opa, c(1,2))
#' getSampleOutlier(bcm.opa, c("GSM1124929","GSM1124941"))

setGeneric("getSampleOutlier",
					 function(profileMatrix, sample.name){ standardGeneric("getSampleOutlier")})


#' @describeIn getSampleOutlier A method for getSampleOutlier
#' @export
setMethod("getSampleOutlier", signature(profileMatrix = "matrix", sample.name = "nuchar"),
					function(profileMatrix, sample.name){
						up.list <- getOutlier(profileMatrix, sample.name)(1)
						down.list <- getOutlier(profileMatrix, sample.name)(-1)



						if(length(sample.name) == 1){
							res.list <- list(up = unlist(up.list, use.names = TRUE), down = unlist(down.list, use.names = TRUE))
						}else{

							res.list <- Map(function(x,y){

								list(up = x, down = y)

							}, up.list, down.list)
							res.list <- unlist(res.list, recursive = FALSE)

						}


          res.list


					})


#' @describeIn getSampleOutlier A method for getSampleOutlier
#' @export
setMethod("getSampleOutlier", signature(profileMatrix = "OPPARList", sample.name = "nuchar"),
					function(profileMatrix, sample.name){
						getSampleOutlier(profileMatrix = profileMatrix@profileMatrix, sample.name = sample.name)
					})









