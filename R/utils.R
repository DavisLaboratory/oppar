# this is a helper function for the below funcntion
get_q <- function(q.list){
	function(x,y){
		q.list[[x]][[y]]
	}

}



getOutlier <- function(profileMatrix, sample.name){
	function(outlr.indicator){
		res <- apply(profileMatrix[, sample.name, drop = FALSE], 2, function(x){ which(x == outlr.indicator)})
		lapply(res,function(x){rownames(profileMatrix[x, , drop = FALSE])})
	}
}





