#' @import GSVA
#' @importFrom Biobase ExpressionSet
#' @import GSEABase
#' @useDynLib oppar matrix_density_R
NULL

##
## function: gsva
## purpose: main function of the package which estimates activity
##          scores for each given gene-set

#' gsva
#'
#' Gene Set Variation Analysis
#' @param expr Gene expression data which can be given either as an \code{ExpressionSet}
#'    object or as a matrix of expression values where rows correspond
#'    to genes and columns correspond to samples.
#' @param gset.idx.list Gene sets provided either as a \code{list} object or as a
#'    \code{GeneSetCollection} object.
#' @param annotation In the case of calling \code{gsva()} with expression data in a \code{matrix}
#'    and gene sets as a \code{GeneSetCollection} object, the \code{annotation} argument
#'    can be used to supply the name of the Bioconductor package that contains
#'    annotations for the class of gene identifiers occurring in the row names of
#'    the expression data matrix. By default \code{gsva()} will try to match the
#'    identifiers in \code{expr} to the identifiers in \code{gset.idx.list} just as
#'    they are, unless the \code{annotation} argument is set.
#' @param method Method to employ in the estimation of gene-set enrichment scores per sample. By default
#'    this is set to \code{gsva} (Hanzelmann et al, 2013) and other options are
#'    \code{ssgsea} (Barbie et al, 2009), \code{zscore} (Lee et al, 2008) or \code{plage}
#'    (Tomfohr et al, 2005). The latter two standardize first expression profiles into z-scores
#'    over the samples and, in the case of \code{zscore}, it combines them together as their sum
#'    divided by the square-root of the size of the gene set,
#'    while in the case of \code{plage} they are used to calculate the singular value decomposition
#'    (SVD) over the genes in the gene set and use the coefficients of the first right-singular vector
#'    as pathway activity profile.
#' @param rnaseq Flag to inform whether the input gene expression data comes from microarray
#'    (\code{rnaseq=FALSE}, default) or RNA-Seq (\code{rnaseq=TRUE}) experiments.
#' @param abs.ranking Flag to determine whether genes should be ranked according to
#'    their sign (\code{abs.ranking=FALSE}) or by absolute value (\code{abs.ranking=TRUE}).
#'    In the latter, pathways with genes enriched on either extreme
#'    (high or low) will be regarded as 'highly' activated.
#' @param min.sz Minimum size of the resulting gene sets.
#' @param max.sz Maximum size of the resulting gene sets.
#' @param no.bootstraps Number of bootstrap iterations to perform.
#' @param bootstrap.percent .632 is the ideal percent samples bootstrapped.
#' @param parallel.sz Number of processors to use when doing the calculations in parallel.
#'    This requires to previously load either the \code{parallel} or the
#'    \code{snow} library. If \code{parallel} is loaded and this argument
#'    is left with its default value (\code{parallel.sz=0}) then it will use
#'    all available core processors unless we set this argument with a
#'    smaller number. If \code{snow} is loaded then we must set this argument
#'    to a positive integer number that specifies the number of processors to
#'    employ in the parallel calculation.
#' @param parallel.type Type of cluster architecture when using \code{snow}.
#' @param mx.diff Offers two approaches to calculate the enrichment statistic (ES)
#'    from the KS random walk statistic. \code{mx.diff=FALSE}: ES is calculated as
#'    the maximum distance of the random walk from 0. \code{mx.diff=TRUE} (default): ES
#'    is calculated as the magnitude difference between the largest positive
#'    and negative random walk deviations.
#' @param tau Exponent defining the weight of the tail in the random walk performed by both the \code{gsva}
#'    (Hanzelmann et al., 2013) and the \code{ssgsea} (Barbie et al., 2009) methods. By default,
#'    this \code{tau=1} when \code{method="gsva"} and \code{tau=0.25} when \code{method="ssgsea"} just
#'    as specified by Barbie et al. (2009) where this parameter is called \code{alpha}.
#' @param kernel Logical, set to \code{TRUE} when the GSVA method employes a kernel non-parametric
#'    estimation of the empirical cumulative distribution function (default) and \code{FALSE}
#'    when this function is directly estimated from the observed data. This last option is
#'    justified in the limit of the size of the sample by the so-called Glivenko-Cantelli theorem.
#' @param ssgsea.norm Logical, set to \code{TRUE} (default) with \code{method="ssgsea"} runs the SSGSEA method
#'    from Barbie et al. (2009) normalizing the scores by the absolute difference between
#'    the minimum and the maximum, as described in their paper. When \code{ssgsea.norm=FALSE}
#'    this last normalization step is skipped.
#' @param verbose Gives information about each calculation step. Default: \code{FALSE}.
#' @param is.gset.list.up.down logical. Is the gene list divided into up/down sublists?
#'
#' @param ... other optional arguments.
#' @export
#' @import GSVA
#' @seealso Hanzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation analysis for
#'     microarray and RNA-Seq data. BMC Bioinformatics, 14, 7. http://doi.org/10.1186/1471-2105-14-7
#' @examples
#' data("Maupin")
#' names(maupin)
#' geneSet<- maupin$sig$EntrezID    #Symbol  ##EntrezID # both up and down genes:
#' up_sig<- maupin$sig[maupin$sig$upDown == "up",]
#' d_sig<- maupin$sig[maupin$sig$upDown == "down",]
#' u_geneSet<- up_sig$EntrezID   #Symbol   # up_sig$Symbol  ## EntrezID
#' d_geneSet<- d_sig$EntrezID
#' #cache(es.dif <- gsva(maupin$data, list(up = u_geneSet, down= d_geneSet), mx.diff=1,
#' #    verbose=TRUE, abs.ranking=FALSE, is.gset.list.up.down=TRUE )$es.obs,
#' #    dir = "data/", prefix = "")
#'
#' ## using cached gsva output
#' data(es.dif)
#' head(es.dif)
setGeneric("gsva", function(expr, gset.idx.list, ...) standardGeneric("gsva"))

#' @describeIn gsva Method for ExpressionSet and list
#' @export
setMethod("gsva", signature(expr="ExpressionSet", gset.idx.list="list"),
          function(expr, gset.idx.list, annotation,
  method=c("gsva", "ssgsea", "zscore", "plage"),
  rnaseq=FALSE,
  abs.ranking=FALSE,
  min.sz=1,
  max.sz=Inf,
  no.bootstraps=0,
  bootstrap.percent = .632,
  parallel.sz=0,
  parallel.type="SOCK",
  mx.diff=TRUE,
  tau=switch(method, gsva=1, ssgsea=0.25, NA),
  kernel=TRUE,
  ssgsea.norm=TRUE,
  verbose=TRUE,
  is.gset.list.up.down = FALSE)
{
  method <- match.arg(method)

  ## filter out genes with constant expression values
  sdGenes <- Biobase::esApply(expr, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            "genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
    }
  }

  if (nrow(expr) < 2)
    stop("Less than two genes in the input ExpressionSet object\n")

  ## map to the actual features for which expression data is available
  mapped.gset.idx.list <- lapply(gset.idx.list,
                                 function(x, y) na.omit(match(x, y)),
                                 Biobase::featureNames(expr))

  if (length(unlist(mapped.gset.idx.list, use.names=FALSE)) == 0)
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")

  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
  mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                         min.sz=max(1, min.sz),
                                         max.sz=max.sz)

  eSco <- .gsva(Biobase::exprs(expr), mapped.gset.idx.list, method, rnaseq, abs.ranking,
                no.bootstraps, bootstrap.percent, parallel.sz, parallel.type,
                mx.diff, tau, kernel, ssgsea.norm, verbose, is.gset.list.up.down )

  if (method != "gsva")
    eSco <- list(es.obs=eSco, bootstrap=NULL, p.vals.sign=NULL)

  eScoEset <- new("ExpressionSet", exprs=eSco$es.obs, phenoData= Biobase::phenoData(expr),
                  experimentData= Biobase::experimentData(expr), annotation="")

	return(list(es.obs=eScoEset,
				      bootstrap=eSco$bootstrap,
              p.vals.sign=eSco$p.vals.sign))
})

#' @describeIn gsva Method for ExpressionSet and GeneSetCollection
#' @export
setMethod("gsva", signature(expr="ExpressionSet", gset.idx.list="GeneSetCollection"),
          function(expr, gset.idx.list, annotation,
  method=c("gsva", "ssgsea", "zscore", "plage"),
  rnaseq=FALSE,
  abs.ranking=FALSE,
  min.sz=1,
  max.sz=Inf,
  no.bootstraps=0,
  bootstrap.percent = .632,
  parallel.sz=0,
  parallel.type="SOCK",
  mx.diff=TRUE,
  tau=switch(method, gsva=1, ssgsea=0.25, NA),
  kernel=TRUE,
  ssgsea.norm=TRUE,
  verbose=TRUE,
  is.gset.list.up.down = FALSE)
{
  method <- match.arg(method)

  ## filter out genes with constant expression values
  sdGenes <- Biobase::esApply(expr, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            "genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
    }
  }

  if (nrow(expr) < 2)
    stop("Less than two genes in the input ExpressionSet object\n")

  if (verbose)
    cat("Mapping identifiers between gene sets and feature names\n")

  ## map gene identifiers of the gene sets to the features in the chip
  mapped.gset.idx.list <- GSEABase::mapIdentifiers(gset.idx.list,
                                                   GSEABase::AnnoOrEntrezIdentifier(Biobase::annotation(expr)))

  ## map to the actual features for which expression data is available
  tmp <- lapply(geneIds(mapped.gset.idx.list),
                                 function(x, y) na.omit(match(x, y)),
                                 Biobase::featureNames(expr))
  names(tmp) <- names(mapped.gset.idx.list)
  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
  mapped.gset.idx.list <- filterGeneSets(tmp,
                                         min.sz=max(1, min.sz),
                                         max.sz=max.sz)

  eSco <- .gsva(Biobase::exprs(expr), mapped.gset.idx.list, method, rnaseq, abs.ranking,
                no.bootstraps, bootstrap.percent, parallel.sz, parallel.type,
                mx.diff, tau, kernel, ssgsea.norm, verbose, is.gset.list.up.down)

  if (method != "gsva")
    eSco <- list(es.obs=eSco, bootstrap=NULL, p.vals.sign=NULL)

  eScoEset <- new("ExpressionSet", exprs=eSco$es.obs, phenoData=Biobase::phenoData(expr),
                  experimentData= Biobase::experimentData(expr), annotation="")

	return(list(es.obs=eScoEset,
				      bootstrap=eSco$bootstrap,
              p.vals.sign=eSco$p.vals.sign))
})

#' @describeIn gsva Method for matrix and GeneSetCollection
#' @export
setMethod("gsva", signature(expr="matrix", gset.idx.list="GeneSetCollection"),
          function(expr, gset.idx.list, annotation,
  method=c("gsva", "ssgsea", "zscore", "plage"),
  rnaseq=FALSE,
  abs.ranking=FALSE,
  min.sz=1,
  max.sz=Inf,
  no.bootstraps=0,
  bootstrap.percent = .632,
  parallel.sz=0,
  parallel.type="SOCK",
  mx.diff=TRUE,
  tau=switch(method, gsva=1, ssgsea=0.25, NA),
  kernel=TRUE,
  ssgsea.norm=TRUE,
  verbose=TRUE,
  is.gset.list.up.down = FALSE)
{
  method <- match.arg(method)

  ## filter out genes with constant expression values
  sdGenes <- apply(expr, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            "genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), , drop=FALSE]
    }
  }

  if (nrow(expr) < 2)
    stop("Less than two genes in the input expression data matrix\n")

  ## map gene identifiers of the gene sets to the features in the matrix
  mapped.gset.idx.list <- gset.idx.list
  if (!missing(annotation)) {
    if (verbose)
      cat("Mapping identifiers between gene sets and feature names\n")

    mapped.gset.idx.list <- GSEABase::mapIdentifiers(gset.idx.list,
                                                     GSEABase::AnnoOrEntrezIdentifier(annotation))
  }

  ## map to the actual features for which expression data is available
  tmp <- lapply(geneIds(mapped.gset.idx.list),
                                 function(x, y) na.omit(match(x, y)),
                                 rownames(expr))
  names(tmp) <- names(mapped.gset.idx.list)

  if (length(unlist(tmp, use.names=FALSE)) == 0)
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")

  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
  mapped.gset.idx.list <- filterGeneSets(tmp,
                                         min.sz=max(1, min.sz),
                                         max.sz=max.sz)

  .gsva(expr, mapped.gset.idx.list, method, rnaseq, abs.ranking,
        no.bootstraps, bootstrap.percent, parallel.sz, parallel.type,
        mx.diff, tau, kernel, ssgsea.norm, verbose, is.gset.list.up.down)
})

#' @describeIn gsva Method for matrix and list
#' @export
setMethod("gsva", signature(expr="matrix", gset.idx.list="list"),
          function(expr, gset.idx.list, annotation,
  method=c("gsva", "ssgsea", "zscore", "plage"),
  rnaseq=FALSE,
  abs.ranking=FALSE,
  min.sz=1,
  max.sz=Inf,
  no.bootstraps=0,
  bootstrap.percent = .632,
  parallel.sz=0,
  parallel.type="SOCK",
  mx.diff=TRUE,
  tau=switch(method, gsva=1, ssgsea=0.25, NA),
  kernel=TRUE,
  ssgsea.norm=TRUE,
  verbose=TRUE,
  is.gset.list.up.down = FALSE)
{
  method <- match.arg(method)

  ## filter out genes with constant expression values
  sdGenes <- apply(expr, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            "genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), , drop=FALSE]
    }
  }

  if (nrow(expr) < 2)
    stop("Less than two genes in the input expression data matrix\n")

  mapped.gset.idx.list <- lapply(gset.idx.list,
                                 function(x ,y) na.omit(match(x, y)),
                                 rownames(expr))

  if (length(unlist(mapped.gset.idx.list, use.names=FALSE)) == 0)
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")

  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
  mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                         min.sz=max(1, min.sz),
                                         max.sz=max.sz)

  .gsva(expr, mapped.gset.idx.list, method, rnaseq, abs.ranking, no.bootstraps,
        bootstrap.percent, parallel.sz, parallel.type,
        mx.diff, tau, kernel, ssgsea.norm, verbose, is.gset.list.up.down )
})

.gsva <- function(expr, gset.idx.list,
  method=c("gsva", "ssgsea", "zscore", "plage"),
  rnaseq=FALSE,
  abs.ranking=FALSE,
  no.bootstraps=0,
  bootstrap.percent = .632,
  parallel.sz=0,
  parallel.type="SOCK",
  mx.diff=TRUE,
  tau=1,
  kernel=TRUE,
  ssgsea.norm=TRUE,
  verbose=TRUE,
  is.gset.list.up.down = FALSE)
{
	if(length(gset.idx.list) == 0){
		stop("The gene set list is empty!  Filter may be too stringent.")
	}

  if (method == "ssgsea") {
	  if(verbose)
		  cat("Estimating ssGSEA scores for", length(gset.idx.list),"gene sets.\n")

    return(ssgsea(expr, gset.idx.list, alpha=tau, parallel.sz=parallel.sz,
                  parallel.type=parallel.type, normalization=ssgsea.norm,
                  verbose=verbose, is.gset.list.up.down))
  }

  if (method == "zscore") {
    if (rnaseq)
      stop("rnaseq=TRUE does not work with method='zscore'.")

	  if(verbose)
		  cat("Estimating combined z-scores for", length(gset.idx.list),"gene sets.\n")

    return(zscore(expr, gset.idx.list, parallel.sz, parallel.type, verbose))
  }

  if (method == "plage") {
    if (rnaseq)
      stop("rnaseq=TRUE does not work with method='plage'.")

	  if(verbose)
		  cat("Estimating PLAGE scores for", length(gset.idx.list),"gene sets.\n")

    return(plage(expr, gset.idx.list, parallel.sz, parallel.type, verbose))
  }

	if(verbose)
		cat("Estimating GSVA scores for", length(gset.idx.list),"gene sets.\n")

	if(parallel.sz > 0 && no.bootstraps > 0){
		if((no.bootstraps %% parallel.sz) != 0){
			stop("'parrallel.sz' must be an integer divisor of 'no.bootsraps'" )
		}
	}
	n.samples <- ncol(expr)
	n.genes <- nrow(expr)
	#n.gset <- length(gset.idx.list)
	n.gset <- ifelse(is.gset.list.up.down, 1, length(gset.idx.list))
	if(is.gset.list.up.down){
		es.obs.row.names <- "GeneSet"
	}else{
		es.obs.row.names <- names(gset.idx.list)
	}


	es.obs <- matrix(NaN, n.gset, n.samples, dimnames=list(es.obs.row.names,colnames(expr)))
	colnames(es.obs) <- colnames(expr)



	rownames(es.obs) <- es.obs.row.names
	#rownames(es.obs) <- names(gset.idx.list)


	if (verbose)
    cat("Computing observed enrichment scores\n")
	es.obs <- compute.geneset.es(expr, gset.idx.list, 1:n.samples,
                               rnaseq=rnaseq, abs.ranking=abs.ranking, parallel.sz=parallel.sz,
                               parallel.type=parallel.type, mx.diff=mx.diff, tau=tau,
                               kernel=kernel, verbose=verbose,
															 is.gset.list.up.down = is.gset.list.up.down)

	# es.bootstraps -> n.gset by n.samples by n.resamples
	es.bootstraps=NULL
	p.vals.wilcoxon=NULL
	p.vals.sign=NULL

	if(no.bootstraps > 0){
		if(verbose) cat("Computing bootstrap enrichment scores\n")
		bootstrap.nsamples <- floor(bootstrap.percent * n.samples)

		p.vals.sign <- matrix(NaN, n.gset, n.samples,
                          dimnames=list(es.obs.row.names, colnames(expr)))

		es.bootstraps <- array(NaN, c(n.gset, n.samples, no.bootstraps))
		if(parallel.sz > 0){

		  if(!.isPackageLoaded("snow")) {
			  stop("Please load the 'snow' library")
		  }
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl <- get("makeCluster", mode="function")
      clSetupRNG <- get("clusterSetupRNG", mode="function")
      clEvalQ <- get("clusterEvalQ", mode="function")
      clExport <- get("clusterExport", mode="function")
      stopCl <- get("stopCluster", mode="function")

			cl <- makeCl(parallel.sz, type = parallel.type)
			.GlobalEnv[["expr"]] <- expr
			.GlobalEnv[["bootstrap.nsamples"]] <- bootstrap.nsamples
			.GlobalEnv[["n.samples"]] <- n.samples
			.GlobalEnv[["gset.idx.list"]] <- gset.idx.list
			clExport(cl,"expr")
			clExport(cl,"bootstrap.nsamples")
			clExport(cl, "n.samples")
			clExport(cl, "gset.idx.list")
			clEvalQ(cl, requireNamespace("oppar", quietly = TRUE))

			clSetupRNG(cl)

			if(verbose) cat("Parallel bootstrap...\n")
			## parallelized bootstrap
			n.cycles <- floor(no.bootstraps / parallel.sz)
			for(i in 1:n.cycles){
				if(verbose) cat("bootstrap cycle ", i, "\n")
				r <- clEvalQ(cl, compute.geneset.es(expr, gset.idx.list,
								sample(n.samples, bootstrap.nsamples, replace=TRUE),
								rnaseq=rnaseq, abs.ranking=abs.ranking, mx.diff=mx.diff,
                tau=tau, kernel=kernel, verbose=verbose, is.gset.list.up.down = is.gset.list.up.down))
				for(j in 1:length(r)){
					es.bootstraps[,,(parallel.sz * (i-1) + j)] <- r[[j]]
				}
			}
			stopCl(cl)
		}else{
			if(verbose) cat("Sequential bootstrap...\n")
			for(i in 1:no.bootstraps){
				es.bootstraps[,,i] <- compute.geneset.es(expr, gset.idx.list,
						sample(n.samples, bootstrap.nsamples, replace=TRUE),
						rnaseq=rnaseq, abs.ranking=abs.ranking, mx.diff=mx.diff,
            tau=tau, kernel=kernel, verbose=verbose, is.gset.list.up.down = is.gset.list.up.down)
			}
		}


		for(i in 1:n.gset){

			for(j in 1:n.samples){
				# non-parametric test if median of empirical dist is 0
				if(es.obs[i,j] > 0){
					p.vals.sign[i,j] <- (1 + sum(es.bootstraps[i,j,] < 0)) / (1 + no.bootstraps)
				}else{
					p.vals.sign[i,j] <- (1 + sum(es.bootstraps[i,j,] > 0)) / (1 + no.bootstraps)
				}
			}
		}
	}

	colnames(es.obs) <- colnames(expr)

	rownames(es.obs) <- es.obs.row.names
	#rownames(es.obs) <- names(gset.idx.list)
	return(list(es.obs=es.obs,
				      bootstrap=list(es.bootstraps=es.bootstraps,
              p.vals.sign=p.vals.sign)))
}


compute.gene.density <- function(expr, sample.idxs, rnaseq=FALSE, kernel=TRUE){
	n.test.samples <- ncol(expr)
	n.genes <- nrow(expr)
	n.density.samples <- length(sample.idxs)

  gene.density <- NA
  if (kernel) {
	  A = .C("matrix_density_R",
			as.double(t(expr[ ,sample.idxs, drop=FALSE])),
			as.double(t(expr)),
			R = double(n.test.samples * n.genes),
			n.density.samples,
			n.test.samples,
			n.genes,
      as.integer(rnaseq))$R

	  gene.density <- t(matrix(A, n.test.samples, n.genes))
  } else {
    gene.density <- t(apply(expr, 1, function(x, sample.idxs) {
                                     f <- ecdf(x[sample.idxs])
                                     f(x)
                                   }, sample.idxs))
    gene.density <- log(gene.density / (1-gene.density))
  }

	return (gene.density)
}

compute.geneset.es <- function(expr, gset.idx.list, sample.idxs, rnaseq=FALSE,
                               abs.ranking, parallel.sz=0, parallel.type="SOCK",
                               mx.diff=TRUE, tau=1, kernel=TRUE, verbose=TRUE,
															 is.gset.list.up.down = FALSE){
	num_genes <- nrow(expr)
	if (verbose) {
    if (kernel) {
      if (rnaseq)
        cat("Estimating ECDFs in rnaseq data with Poisson kernels\n")
      else
        cat("Estimating ECDFs in microarray data with Gaussian kernels\n")
    } else
      cat("Estimating ECDFs directly\n")
  }
	gene.density <- compute.gene.density(expr, sample.idxs, rnaseq, kernel)

	compute_rank_score <- function(sort_idx_vec){
		tmp <- rep(0, num_genes)
		tmp[sort_idx_vec] <- abs(seq(from=num_genes,to=1) - num_genes/2)
		return (tmp)
	}

	rank.scores <- rep(0, num_genes)
	if(abs.ranking){
		sort.sgn.idxs <- apply(abs(gene.density), 2, order, decreasing=TRUE)# n.genes * n.samples
		if(is.gset.list.up.down){
			sort.sgn.idxs.down <-  apply(abs(gene.density), 2, order, decreasing=FALSE)}
	}else{
		sort.sgn.idxs <- apply(gene.density, 2, order, decreasing=TRUE) # n.genes * n.samples
		if(is.gset.list.up.down){
			sort.sgn.idxs.down <- apply(gene.density, 2, order, decreasing=FALSE)
		}
	}

	rank.scores <- apply(sort.sgn.idxs, 2, compute_rank_score)
	if(is.gset.list.up.down){
		rank.scores.divided.gset.list <- apply(sort.sgn.idxs.down, 2, compute_rank_score)
	}
	haveParallel <- .isPackageLoaded("parallel")
	haveSnow <- .isPackageLoaded("snow")

	if (parallel.sz > 1 || haveParallel) {
		if (!haveParallel && !haveSnow) {
			stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
		}

    if (!haveParallel) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl <- get("makeCluster", mode="function")
      parSapp <- get("parSapply", mode="function")
      clEvalQ <- get("clusterEvalQ", mode="function")
      stopCl <- get("stopCluster", mode="function")

      if (verbose)
        cat("Allocating cluster\n")
		  cl <- makeCl(parallel.sz, type = parallel.type)
		  clEvalQ(cl, requireNamespace("oppar", quietly = TRUE))
		  if (verbose) {
		    cat("Estimating enrichment scores in parallel\n")
	      if(mx.diff) {
          cat("Taking diff of max KS.\n")
        } else{
          cat("Evaluting max KS.\n")
        }
		  }

		  if(is.gset.list.up.down){
		  	# can we perhaps re-write this with Reduce
		  	m <- t(parSapp(cl, gset.idx.list["up"], ks_test_m,
		  								 gene.density=rank.scores,
		  								 sort.idxs=sort.sgn.idxs,
		  								 mx.diff=mx.diff, tau=tau, verbose=FALSE))

		  	m.down <- t(parSapp(cl, gset.idx.list["down"], ks_test_m,
		  								 gene.density=rank.scores.divided.gset.list,
		  								 sort.idxs=sort.sgn.idxs.down,
		  								 mx.diff=mx.diff, tau=tau, verbose=FALSE))
		  	m <- m + m.down
		  }else{
		  	m <- t(parSapp(cl, gset.idx.list, ks_test_m,
		  								 gene.density=rank.scores,
		  								 sort.idxs=sort.sgn.idxs,
		  								 mx.diff=mx.diff, tau=tau, verbose=FALSE))
		  	}




		  if(verbose)
        cat("Cleaning up\n")
		  stopCl(cl)

    } else {             ## use parallel

      mclapp <- get('mclapply', envir=getNamespace('parallel'))
      detCor <- get('detectCores', envir=getNamespace('parallel'))
      nCores <- detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)

      pb <- NULL
      if (verbose){
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
        assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
        assign("nGeneSets", ceiling(length(gset.idx.list) / getOption("mc.cores")), envir=globalenv())
        assign("iGeneSet", 0, envir=globalenv())
      }


      if(is.gset.list.up.down){
      	m <- mclapp(gset.idx.list["up"], ks_test_m,
      							gene.density=rank.scores,
      							sort.idxs=sort.sgn.idxs,
      							mx.diff=mx.diff, tau=tau, verbose=verbose)
      	m <- do.call("rbind", m)
      	colnames(m) <- colnames(expr)

      	m.down <- mclapp(gset.idx.list["down"], ks_test_m,
      									 gene.density=rank.scores.divided.gset.list,
      									 sort.idxs=sort.sgn.idxs.down,
      									 mx.diff=mx.diff, tau=tau, verbose=verbose)
      	m.down <- do.call("rbind", m.down)
      	colnames(m.down) <- colnames(expr)

      	m <- m + m.down

      }else{
      	m <- mclapp(gset.idx.list, ks_test_m,
      							gene.density=rank.scores,
      							sort.idxs=sort.sgn.idxs,
      							mx.diff=mx.diff, tau=tau, verbose=verbose)
      	m <- do.call("rbind", m)
      	colnames(m) <- colnames(expr)
      }



      if (verbose) {
        close(get("progressBar", envir=globalenv()))
      }
    }

	} else {
		if (verbose) {
      cat("Estimating enrichment scores\n")
	    if (mx.diff) {
        cat("Taking diff of max KS.\n")
      } else{
        cat("Evaluting max KS.\n")
      }
    }
    pb <- NULL
    if (verbose){
      assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
      assign("nGeneSets", length(gset.idx.list), envir=globalenv())
      assign("iGeneSet", 0, envir=globalenv())
    }


		if(is.gset.list.up.down){
			m <- t(sapply(gset.idx.list["up"], ks_test_m, rank.scores, sort.sgn.idxs,
										mx.diff=mx.diff, tau=tau, verbose=verbose))
		  m.down <- t(sapply(gset.idx.list["down"], ks_test_m, rank.scores.divided.gset.list,
											 sort.sgn.idxs.down,
											 mx.diff=mx.diff, tau=tau, verbose=verbose))

		  m <- m + m.down
		}else{
			m <- t(sapply(gset.idx.list, ks_test_m, rank.scores, sort.sgn.idxs,
										mx.diff=mx.diff, tau=tau, verbose=verbose))
		}


    if (verbose) {
      setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
      close(get("progressBar", envir=globalenv()))
    }
	}
	return (m)
}


ks_test_m <- function(gset_idxs, gene.density, sort.idxs, mx.diff=TRUE, tau=1, verbose=TRUE){

	n.genes <- nrow(gene.density)
	n.samples <- ncol(gene.density)
	n.geneset <- length(gset_idxs)

	geneset.sample.es = .C("ks_matrix_R",
			as.double(gene.density),
			R = double(n.samples),
			as.integer(sort.idxs),
			n.genes,
			as.integer(gset_idxs),
			n.geneset,
			as.double(tau),
			n.samples,
			as.integer(mx.diff))$R

  if (verbose) {
    assign("iGeneSet", get("iGeneSet", envir=globalenv()) + 1, envir=globalenv())
    setTxtProgressBar(get("progressBar", envir=globalenv()),
                      get("iGeneSet", envir=globalenv()) / get("nGeneSets", envir=globalenv()))
  }

	return (geneset.sample.es)
}


## ks-test in R code - testing only
ks_test_Rcode <- function(gene.density, gset_idxs, tau=1, make.plot=FALSE){

	n.genes = length(gene.density)
	n.gset = length(gset_idxs)

	sum.gset <- sum(abs(gene.density[gset_idxs])^tau)

	dec = 1 / (n.genes - n.gset)

	sort.idxs <- order(gene.density,decreasing=TRUE)
	offsets <- sort(match(gset_idxs, sort.idxs))

	last.idx = 0
	values <- rep(NaN, length(gset_idxs))
	current = 0
	for(i in seq_along(offsets)){
		current = current + abs(gene.density[sort.idxs[offsets[i]]])^tau / sum.gset - dec * (offsets[i]-last.idx-1)

		values[i] = current
		last.idx = offsets[i]
	}
	check_zero = current - dec * (n.genes-last.idx)
	#if(check_zero > 10^-15){
	#	stop(paste=c("Expected zero sum for ks:", check_zero))
	#}
	if(make.plot){ plot(offsets, values,type="l") }

	max.idx = order(abs(values),decreasing=TRUE)[1]
	mx.value <- values[max.idx]

	return (mx.value)
}

rndWalk <- function(gSetIdx, geneRanking, j, R, alpha) {
  indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
  stepCDFinGeneSet <- cumsum((abs(R[geneRanking, j]) *
                      indicatorFunInsideGeneSet)^alpha) /
                      sum((abs(R[geneRanking, j]) *
                      indicatorFunInsideGeneSet)^alpha)
  stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet) /
                       sum(!indicatorFunInsideGeneSet)
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

  sum(walkStat)
}

ssgsea <- function(X, geneSets, alpha=0.25, parallel.sz,
                   parallel.type, normalization=TRUE, verbose, is.gset.list.up.down) {

  p <- nrow(X)
  n <- ncol(X)

  if (verbose) {
    assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
    assign("nSamples", n, envir=globalenv())
    assign("iSample", 0, envir=globalenv())
  }

  R <- apply(X, 2, function(x,p) as.integer(rank(x)), p)

	haveParallel <- .isPackageLoaded("parallel")
	haveSnow <- .isPackageLoaded("snow")

  cl <- makeCl <- parSapp <- stopCl <- mclapp <- detCor <- nCores <- NA
	if (parallel.sz > 1 || haveParallel) {
		if (!haveParallel && !haveSnow) {
			stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
		}

    if (!haveParallel) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl <- get("makeCluster", mode="function")
      parSapp <- get("parSapply", mode="function")
      stopCl <- get("stopCluster", mode="function")

      if (verbose)
        cat("Allocating cluster\n")
		  cl <- makeCl(parallel.sz, type = parallel.type)
    } else {             ## use parallel

      mclapp <- get('mclapply', envir=getNamespace('parallel'))
      detCor <- get('detectCores', envir=getNamespace('parallel'))
      nCores <- detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)
      if (verbose)
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
    }
  }

  es <- sapply(1:n, function(j, R, geneSets, alpha, is.gset.list.up.down ) {
                      if (verbose) {
                        assign("iSample", get("iSample", envir=globalenv()) + 1, envir=globalenv())
                        setTxtProgressBar(get("progressBar", envir=globalenv()),
                                          get("iSample", envir=globalenv()) / get("nSamples", envir=globalenv()))
                      }

                      geneRanking <- order(R[, j], decreasing=TRUE)
                      geneRanking.down <- order(R[, j], decreasing=FALSE)
                      es_sample <- NA
                      if (parallel.sz == 1 || (is.na(cl) && !haveParallel)){
                      	if(is.gset.list.up.down){
                      		es_sample_up <- unlist(sapply(geneSets["up"], rndWalk, geneRanking, j, R, alpha))
                      		es_sample_down <- unlist(sapply(geneSets["down"], rndWalk, geneRanking.down, j, R, alpha))
                      		es_sample <- es_sample_up + es_sample_down


                      	}else{
                        es_sample <- unlist(sapply(geneSets, rndWalk, geneRanking, j, R, alpha))}

                      }else {

                        if (is.na(cl)){
                        	if(is.gset.list.up.down){
                        		es_sample_up <- unlist(mclapp(geneSets["up"], rndWalk, geneRanking, j, R, alpha))
                        		es_sample_down <- unlist(mclapp(geneSets["down"], rndWalk, geneRanking.down, j, R, alpha))
                        		es_sample <- es_sample_up + es_sample_down

                        	}else{
                          es_sample <- unlist(mclapp(geneSets, rndWalk, geneRanking, j, R, alpha))}

                        }else{
                        	if(is.gset.list.up.down){
                        		print("is.working")
                        		es_sample_up <- unlist(parSapp(cl, geneSets["up"], rndWalk, geneRanking, j, R, alpha))
                        		es_sample_down <- unlist(parSapp(cl, geneSets["down"], rndWalk, geneRanking.down, j, R, alpha))
                        		es_sample <- es_sample_up + es_sample_down

                        	}else{
                          es_sample <- unlist(parSapp(cl, geneSets, rndWalk, geneRanking, j, R, alpha))}
                        }

                     }

                      es_sample
                    }, R, geneSets, alpha, is.gset.list.up.down)

# fix rownames

  if (length(geneSets) == 1 || is.gset.list.up.down)
    es <- matrix(es, nrow=1)

  if (normalization) {
    ## normalize enrichment scores by using the entire data set, as indicated
    ## by Barbie et al., 2009, online methods, pg. 2
    es <- apply(es, 2, function(x, es) x / (range(es)[2] - range(es)[1]), es)
  }

  if (length(geneSets) == 1 || is.gset.list.up.down)
    es <- matrix(es, nrow=1)

  rownames(es) <- ifelse(is.gset.list.up.down, "GeneSet", names(geneSets))
  colnames(es) <- colnames(X)

  if (verbose) {
    setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
    close(get("progressBar", envir=globalenv()))
  }

  if (!is.na(cl))
    stopCl(cl)

  es
}

combinez <- function(gSetIdx, j, Z) sum(Z[gSetIdx, j]) / sqrt(length(gSetIdx))

zscore <- function(X, geneSets, parallel.sz, parallel.type, verbose) {

  p <- nrow(X)
  n <- ncol(X)

  if (verbose) {
    assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
    assign("nSamples", n, envir=globalenv())
    assign("iSample", 0, envir=globalenv())
  }

  Z <- t(apply(X, 1, function(x) (x-mean(x))/sd(x)))

	haveParallel <- .isPackageLoaded("parallel")
	haveSnow <- .isPackageLoaded("snow")

  cl <- makeCl <- parSapp <- stopCl <- mclapp <- detCor <- nCores <- NA
	if (parallel.sz > 1 || haveParallel) {
		if (!haveParallel && !haveSnow) {
			stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
		}

    if (!haveParallel) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl <- get("makeCluster", mode="function")
      parSapp <- get("parSapply", mode="function")
      stopCl <- get("stopCluster", mode="function")

      if (verbose)
        cat("Allocating cluster\n")
		  cl <- makeCl(parallel.sz, type = parallel.type)
    } else {             ## use parallel

      mclapp <- get('mclapply', envir=getNamespace('parallel'))
      detCor <- get('detectCores', envir=getNamespace('parallel'))
      nCores <- detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)
      if (verbose)
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
    }
  }

  es <- sapply(1:n, function(j, Z, geneSets) {
                      if (verbose) {
                        assign("iSample", get("iSample", envir=globalenv()) + 1, envir=globalenv())
                        setTxtProgressBar(get("progressBar", envir=globalenv()),
                                          get("iSample", envir=globalenv()) / get("nSamples", envir=globalenv()))
                      }
                      es_sample <- NA
                      if (parallel.sz == 1 || (is.na(cl) && !haveParallel))
                        es_sample <- sapply(geneSets, combinez, j, Z)
                      else {
                        if (is.na(cl))
                          es_sample <- mclapp(geneSets, combinez, j, Z)
                        else
                          es_sample <- parSapp(cl, geneSets, combinez, j, Z)
                      }

                      unlist(es_sample)
                    }, Z, geneSets)

  if (length(geneSets) == 1)
    es <- matrix(es, nrow=1)

  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)

  if (verbose) {
    setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
    close(get("progressBar", envir=globalenv()))
  }

  if (!is.na(cl))
    stopCl(cl)

  es
}

rightsingularsvdvectorgset <- function(gSetIdx, Z) {
    s <- svd(Z[gSetIdx, ])
  s$v[, 1]
}

plage <- function(X, geneSets, parallel.sz, parallel.type, verbose) {

  p <- nrow(X)
  n <- ncol(X)

  if (verbose) {
    assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
    assign("nGeneSets", length(geneSets), envir=globalenv())
    assign("iGeneSet", 0, envir=globalenv())
  }

  Z <- t(apply(X, 1, function(x) (x-mean(x))/sd(x)))

	haveParallel <- .isPackageLoaded("parallel")
	haveSnow <- .isPackageLoaded("snow")

  ## the masterDescriptor() calls are disabled since they are not available in windows
  ## they would help to report progress by just one of the processors. now all processors
  ## will reporting progress. while this might not be the right way to report progress in
  ## parallel it should not affect a correct execution and progress should be more or less
  ## being reported to some extent.
  cl <- makeCl <- parSapp <- stopCl <- mclapp <- detCor <- nCores <- NA ## masterDesc <- NA
	if(parallel.sz > 1 || haveParallel) {
		if(!haveParallel && !haveSnow) {
			stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
		}

    if (!haveParallel) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl <- get("makeCluster", mode="function")
      parSapp <- get("parSapply", mode="function")
      stopCl <- get("stopCluster", mode="function")

      if (verbose)
        cat("Allocating cluster\n")
		  cl <- makeCl(parallel.sz, type = parallel.type)
    } else {             ## use parallel

      mclapp <- get('mclapply', envir=getNamespace('parallel'))
      detCor <- get('detectCores', envir=getNamespace('parallel'))
      ## masterDesc <- get('masterDescriptor', envir=getNamespace('parallel'))
      nCores <- detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)
      if (verbose)
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
    }
  }

  if (parallel.sz == 1 || (is.na(cl) && !haveParallel))
    es <- t(sapply(geneSets, function(gset, Z) {
                             if (verbose) {
                               assign("iGeneSet", get("iGeneSet", envir=globalenv()) + 1, envir=globalenv())
                               setTxtProgressBar(get("progressBar", envir=globalenv()),
                                                 get("iGeneSet", envir=globalenv()) / get("nGeneSets", envir=globalenv()))
                             }
                             rightsingularsvdvectorgset(gset, Z)
                           }, Z))
  else {
    if (is.na(cl)) {
      ## firstproc <- mclapp(as.list(1:(options("mc.cores")$mc.cores)), function(x) masterDesc())[[1]]
      es <- mclapp(geneSets, function(gset, Z) { ##, firstproc) {
                                 if (verbose) { ## && masterDesc() == firstproc) {
                                   assign("iGeneSet", get("iGeneSet", envir=globalenv()) + 1, envir=globalenv())
                                   setTxtProgressBar(get("progressBar", envir=globalenv()),
                                                     get("iGeneSet", envir=globalenv()) / get("nGeneSets", envir=globalenv()))
                                 }
                                 rightsingularsvdvectorgset(gset, Z)
                               }, Z) ##, firstproc)
      es <- do.call(rbind, es)
    } else {
      if (verbose)
        message("Progress reporting for plage with a snow cluster not yet implemented")

      es <- parSapp(geneSets, function(gset, Z) {
                                  if (verbose) {
                                    assign("iGeneSet", get("iGeneSet", envir=globalenv()) + 1, envir=globalenv())
                                    setTxtProgressBar(get("progressBar", envir=globalenv()),
                                                      get("iGeneSet", envir=globalenv()) / get("nGeneSets", envir=globalenv()))
                                  }
                                  rightsingularsvdvectorgset(gset, Z)
                                }, Z)
      es <- do.call(rbind, es)
    }
  }

  if (length(geneSets) == 1)
    es <- matrix(es, nrow=1)

  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)

  if (verbose) {
    setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
    close(get("progressBar", envir=globalenv()))
  }

  if (!is.na(cl))
    stopCl(cl)

  es
}

setGeneric("filterGeneSets", function(gSets, ...) standardGeneric("filterGeneSets"))

setMethod("filterGeneSets", signature(gSets="list"),
          function(gSets, min.sz=1, max.sz=Inf) {
	gSetsLen <- sapply(gSets,length)
	return (gSets[gSetsLen >= min.sz & gSetsLen <= max.sz])
})

setMethod("filterGeneSets", signature(gSets="GeneSetCollection"),
          function(gSets, min.sz=1, max.sz=Inf) {
	gSetsLen <- sapply(geneIds(gSets),length)
	return (gSets[gSetsLen >= min.sz & gSetsLen <= max.sz])
})



setGeneric("computeGeneSetsOverlap", function(gSets, uniqGenes=unique(unlist(gSets, use.names=FALSE)), ...) standardGeneric("computeGeneSetsOverlap"))

setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="character"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  totalGenes <- length(uniqGenes)

  ## map to the features requested
  gSets <- lapply(gSets, function(x, y) as.vector(na.omit(match(x, y))), uniqGenes)

  lenGsets <- sapply(gSets, length)
  totalGsets <- length(gSets)

  gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members <- cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] <- 1

  .computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})

setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  uniqGenes <- Biobase::featureNames(uniqGenes)
  totalGenes <- length(uniqGenes)

  ## map to the actual features for which expression data is available
  gSets <- lapply(gSets, function(x, y) as.vector(na.omit(match(x, y))), uniqGenes)

  lenGsets <- sapply(gSets, length)
  totalGsets <- length(gSets)

  gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members <- cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] <- 1

  .computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})

setMethod("computeGeneSetsOverlap", signature(gSets="GeneSetCollection", uniqGenes="character"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  .computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})

setMethod("computeGeneSetsOverlap", signature(gSets="GeneSetCollection", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  ## map gene identifiers of the gene sets to the features in the chip
  gSets <- GSEABase::mapIdentifiers(gSets, GSEABase::AnnoOrEntrezIdentifier(Biobase::annotation(uniqGenes)))

  uniqGenes <- Biobase::featureNames(uniqGenes)

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  .computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})

.computeGeneSetsOverlap <- function(gSetsMembershipMatrix, min.sz=1, max.sz=Inf) {
  ## gSetsMembershipMatrix should be a (genes x gene-sets) incidence matrix

  lenGsets <- colSums(gSetsMembershipMatrix)

  szFilterMask <- lenGsets >= max(1, min.sz) & lenGsets <= max.sz
  if (!any(szFilterMask))
    stop("No gene set meets the minimum and maximum size filter\n")

  gSetsMembershipMatrix <- gSetsMembershipMatrix[, szFilterMask]
  lenGsets <- lenGsets[szFilterMask]

  totalGsets <- ncol(gSetsMembershipMatrix)

  M <- t(gSetsMembershipMatrix) %*% gSetsMembershipMatrix

  M1 <- matrix(lenGsets, nrow=totalGsets, ncol=totalGsets,
               dimnames=list(colnames(gSetsMembershipMatrix), colnames(gSetsMembershipMatrix)))
  M2 <- t(M1)
  M.min <- matrix(0, nrow=totalGsets, ncol=totalGsets)
  M.min[M1 < M2] <- M1[M1 < M2]
  M.min[M2 <= M1] <- M2[M2 <= M1]
  overlapMatrix <- M / M.min

  return (overlapMatrix)
}

## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: isPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise

.isPackageLoaded <- function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
  (paste("package:", name, sep="") %in% search()) ||
  (name %in% loadedNamespaces())
}

##
## ARE THESE FUNCTIONS STILL NECESSARY ?????
##

##a <- replicate(1000, compute.null.enrichment(10000,50,make.plot=F))

compute.null.enrichment <- function(n.genes, n.geneset, make.plot=FALSE){
	ranks <- (n.genes/2) - rev(1:n.genes)
	#null.gset.idxs <- seq(1, n.genes, by=round(n.genes / n.geneset))
	null.gset.idxs <- sample(n.genes, n.geneset)
	null.es <- ks_test_Rcode(ranks, null.gset.idxs,make.plot=make.plot)
	return (null.es)
}


load.gmt.data <- function(gmt.file.path){
	tmp <- readLines(gmt.file.path)
	gsets <- list()
	for(i in 1:length(tmp)){
		t <- strsplit(tmp[i],'\t')[[1]]
		gsets[[t[1]]] <- t[3:length(t)]
	}
	return (gsets)
}

compute.gset.overlap.score <- function(gset.idxs){
	n <- length(gset.idxs)
	mx.idx <- max(unlist(gset.idxs, use.names=FALSE))
	l <- c(sapply(gset.idxs, length))

	gset.M <- matrix(0, nrow=mx.idx, ncol=n)
	for(i in 1:n){
		gset.M[gset.idxs[[i]],i] = 1
	}
	M <- t(gset.M) %*% gset.M

	M1 <- matrix(l, nrow=n, ncol=n)
	M2 <- t(M1)
	M.min <- matrix(0, nrow=n, ncol=n)
	M.min[M1 < M2] <- M1[M1 < M2]
	M.min[M2 <= M1] <- M2[M2 <= M1]
	M.score <- M / M.min
	return (M.score)
}
