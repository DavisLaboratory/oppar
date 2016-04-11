
# mCOPA generic and method


#' Outlier profile Analysis
#'
#' Returns a matrix with 0, -1 and 1 entries that describe outlier profiles
#' in samples. The rows reperesent genes and the columns represent samples. -1 implies
#' that the gene is a down-regulated outlier, 1 indicates an up-regulate outlier and 0
#'  means that the gene is not an outlier in a sample.
#' @param exprs.matrix Gene expression data. Can be either a matrix or an
#'   object of type ExpressionSet.
#' @param group A vector of factors representing the groups to which each sample belong.
#'   This can be either a vector of 0s and 1s, or normal and cases.
#' @param lower.quantile Numeric. The cut-off for lower quantile when determining outliers.
#'   Default to 0.05
#' @param  upper.quantile Numeric. The cut-off for upper quantile when determining outliers.
#'   Default to 0.95
#' @param ... Numeric. To supply values for upper.quantile and lower.quantile
#'   arguments if default values are going to be override.
#' @return \code{opa} returns an object of type \code{OPPARList}. The outlier profiles
#'   are stored in \code{profileMatrix} and can be accessed using $. It it also
#'   possible to retrieve parameters used to run the outlier profile analysis, such
#'   as \code{upper.quantile}, \code{lower.quantile} via the $ operator.
#' @export
#' @examples
#' # loading bcm object from GSE46141 dataset
#' data(GSE46141)
#' library(Biobase)
#' # defining the group variable. local breast tumors are the controls
#' # and the rest of the samples are the diseased samples
#' group <- sapply(pData(bcm)$source_name_ch1, function(x){ ifelse(x == "breast",0,1)})
#' group <- factor(group)
#' # running opa with default values (i.e upper.quantile = 0.95, lower.quantile = 0.05)
#' # the result is an object of type OPPARList
#' opa(bcm,group = group)
#' @seealso Wang, C., Taciroglu, A., Maetschke, S. R., Nelson, C. C., Ragan, M. A., & Davis, M. J. (2012).
#'     mCOPA: analysis of heterogeneous features in cancer expression data. Journal
#'     of Clinical Bioinformatics, 2, 22. http://doi.org/10.1186/2043-9113-2-22
setGeneric("opa",
					 function(exprs.matrix, ...){
					 	standardGeneric("opa")})


#' @describeIn opa opa(exprs.matrix, group, lower.quantile = 0.05, upper.quantile = 0.95)
#' @export
setMethod("opa", signature(exprs.matrix = "matrix"),
					function(exprs.matrix,group ,upper.quantile = 0.95, lower.quantile = 0.05){
						if (length(group) != ncol(exprs.matrix)) {
							stop(sprintf("number of elements in group doesn't match the number of columns. Expected %d, received %d",
													 ncol(exprs.matrix),length(group)))}
						if (!all(levels(group) == c("0","1"))) levels(group) = c("0","1")

	outlier_profile_matrix <- apply(exprs.matrix, 1, FUN=function(x){
		MAD = mad(x,na.rm = TRUE)
		#MAD = median(abs(x-median(x)))
		if (MAD == 0) {
			exprs.transformed <- rep(0,length(x))
		}else{
			exprs.transformed <- (x-median(x))/MAD

			quantile.vals <- vapply(split(exprs.transformed,group),
															quantile, probs = c(seq(0.25, 1, 0.25),
																									upper.quantile, lower.quantile),
															na.rm = TRUE, FUN.VALUE = numeric(6))

			up.outlier.thr <- quantile.vals["75%","1"] + 1.5 * IQR(exprs.transformed[which(group == 1)], na.rm = TRUE)
			down.outlier.thr <- quantile.vals["25%","1"] - 1.5 * IQR(exprs.transformed[which(group == 1)], na.rm = TRUE)

			# adjustment value if any of the quantiles is 0
			adj.val <- ifelse(any(as.vector(quantile.vals[5:6,]) == 0),  0.1, 0)
			down.fc <- log2(abs((quantile.vals[paste0(lower.quantile * 100, "%"),"1"] + adj.val)/(quantile.vals[paste0(lower.quantile * 100, "%"), "0"] + adj.val)))
			up.fc <- log2(abs((quantile.vals[paste0(upper.quantile * 100, "%"), "1"] + adj.val)/(quantile.vals[paste0(upper.quantile * 100, "%"), "0"] + adj.val)))

			# --------- There is a problem here ------------ #

			# up outliers -------------------------------------
			is.up.outlier.cspls <- exprs.transformed[which(group == 1)] > up.outlier.thr
			null.vec.up.cspls <- logical(length(group))
			null.vec.up.cspls[which(group == 1)] <- is.up.outlier.cspls
			is.up.outlier.cspls <- null.vec.up.cspls

			is.up.outlier.nspls <- exprs.transformed[which(group == 0)] > up.outlier.thr
			null.vec.up.nspls <- logical(length(group))
			null.vec.up.nspls[which(group == 0)] <- is.up.outlier.nspls
			is.up.outlier.nspls <- null.vec.up.nspls

			is.sig.up <- abs(up.fc) >= 2

			# down outliers --------------------------------
			is.down.outlier.cspls <- exprs.transformed[which(group == 1)] < down.outlier.thr
			null.vec.down.cspls <- logical(length(group))
			null.vec.down.cspls[which(group == 1)] <- is.down.outlier.cspls
			is.down.outlier.cspls <- null.vec.down.cspls

			is.down.outlier.nspls <- exprs.transformed[which(group == 0)] < down.outlier.thr
			null.vec.down.nspls <- logical(length(group))
			null.vec.down.nspls[which(group == 0)] <- is.down.outlier.nspls
			is.down.outlier.nspls <- null.vec.down.nspls

			is.sig.down <- abs(down.fc) >= 2
			up.filtered.index <- lapply(list(is.up.outlier.cspls, !is.up.outlier.nspls), which)

			if(all(((sum(is.up.outlier.nspls, na.rm = TRUE) + sum(is.up.outlier.cspls, na.rm = TRUE)) >= 1), is.sig.up )){
				up.filter.index <- Reduce(intersect,up.filtered.index)
				exprs.transformed[up.filter.index] <- 1

			}

			if(all(((sum(is.down.outlier.nspls, na.rm = TRUE) + sum(is.down.outlier.cspls, na.rm = TRUE)) >= 1),is.sig.down)){
				down.filtered.index <- sapply(list(is.down.outlier.cspls, !is.down.outlier.nspls), which )
				exprs.transformed[Reduce(intersect, down.filtered.index)] <- -1


			}


			exprs.transformed[ exprs.transformed != 1 & exprs.transformed != -1] <- 0

		}

		exprs.transformed


	})

	zero_features <- apply(outlier_profile_matrix, 2, function(x){
		ifelse(all(x == 0), TRUE, FALSE)
	})

	outlier_profile_matrix <- outlier_profile_matrix[, -which(zero_features)]
	message(sprintf("%d features had zero entries and were removed from the profile matrix",
									sum(zero_features)))

	outlier_profile_matrix <- outlier_profile_matrix[which(group == 1), ] # change here  to 1 later
	new("OPPARList", profileMatrix = t(outlier_profile_matrix),
			upper.quantile = upper.quantile, lower.quantile = lower.quantile,
			group = group[which(group == 1)])

 })


# method for signature exprs = ExpressionSet


#' @describeIn opa opa(eset, group, lower.quantile = 0.05, upper.quantile = 0.95)
#' @export
#' @importClassesFrom Biobase ExpressionSet
setMethod("opa", signature(exprs.matrix = "ExpressionSet"),
					function(exprs.matrix, group, upper.quantile = 0.95, lower.quantile = 0.05){
						opa(exprs.matrix = Biobase::exprs(exprs.matrix), group, upper.quantile, lower.quantile)
					})












