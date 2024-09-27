
# random row of matrix
multi_sample_permutation <- function(matrix, seed = 1234) {
  set.seed(seed)
  perm_idx <- sample(seq_len(nrow(expression)), size = nrow(expression), replace = FALSE)
  perm_expression <- expression[perm_idx, ]
  rownames(perm_expression) <- colnames(perm_expression) <- NULL
  return(perm_expression)
}

# superimposed matrix
superpositionRank <- function(matrix){
  m <- nrow(matrix)
  n <- ncol(matrix)
  max_val <- max(matrix)

  max_matrix <- matrix(rep(seq(0, max_val*(n-1), by = max_val), each=m), nrow=m)
  superposition <- matrix + max_matrix

  superposition_rank <- matrix(rank(superposition, ties.method = "first"), nrow = m)

  return(superposition_rank)
}

# rankMatrixByCol
rankMatrixByCol <- function(superposition_rank){
  m <- nrow(superposition_rank)
  n <- ncol(superposition_rank)

  minus_matrix <- matrix(rep(seq(0, m*(n-1), by=m), each=m), nrow=m)
  subtract_matrix <- superposition_rank - minus_matrix

  return(subtract_matrix)
}


#' sort matrix by column
#'
#' @description Sort every column of a matrix.
#'
#' @param  raw_matrix A matrix consist of numeric value.
#' @param  decreasing Logical, if true sort is decreasing, default by TRUE.
#'
#'

sortMatrixByCol <- function(raw_matrix, decreasing=TRUE){
  if(is.null(dim(raw_matrix))) {
    subtract_matrix <- matrix(raw_matrix, nrow = 1)
  }else{
    n.num <- ncol(raw_matrix)
    m.num <- nrow(raw_matrix)
    #superposition
    over_matrix <- matrix(rep(seq(0, max(raw_matrix)*(n.num-1), by = max(raw_matrix)), each=m.num), nrow=m.num)

    if(decreasing){
      superposition <- raw_matrix - over_matrix
      superposition_sort <- sort(superposition, decreasing=decreasing)
      subtract_matrix <- superposition_sort + over_matrix
    }else{
      superposition <- raw_matrix + over_matrix
      superposition_sort <- sort(superposition, decreasing=decreasing)
      subtract_matrix <- superposition_sort - over_matrix
    }
  }
  return(subtract_matrix)
}



#' identical gene position
#' @description With the rank of all values in the matrix(superposition.rank),
#' get the matrix consist with 0 and 1, indicated the position of geneset, required for calculating ES score of the geneset.
#'
#' @param  superposition_rank A matrix indicating rank of all values in the raw matrix.
#' @param  gSet_pos_matrix A matrix consist of '0' and '1' show the position of signature genes in expression profile,
#' with same size of expression profile, '1' stands for signature genes expression, '0' for other genes.
#'
#' @return matrix
#'
indicatorInsideGeneSetMatrix <- function(superposition_rank, gSet_pos_matrix){
  m <- nrow(superposition_rank)
  n <- ncol(superposition_rank)

  indicatorFunInsideGeneSet <- numeric(length = m*n)
  indicatorFunInsideGeneSet[superposition_rank] <- gSet_pos_matrix

  indicatorFunInsideGeneSet <- matrix(indicatorFunInsideGeneSet, nrow = m)

  indicatorFunInsideGeneSet <- apply(indicatorFunInsideGeneSet, 2, rev)

  return(indicatorFunInsideGeneSet)

}


#' ES score
#' @description Calculating ES score of one geneset.
#' @param  exp_rank A matrix indicating rank of every element in column.
#' @param  superposition_rank A matrix indicating rank of all values in the matrix.
#' @param  gSet_pos_matrix A matrix consist of '0' and '1' show the position of signature genes in expression profile,
#' with same size of expression profile, '1' stands for signature genes expression, '0' for other genes.
#' @param  alpha Quantity alpha usually set to 0.25.
#'
#' @return vector of numeric values indicating ES score of one geneset for every sample
#'
ssGSEA_ES <- function(exp_rank, superposition_rank, gSet_pos_matrix, alpha=0.25) {

  indicatorInsideGeneSet_matrix <- indicatorInsideGeneSetMatrix(superposition_rank = superposition_rank,
                                                                gSet_pos_matrix = gSet_pos_matrix)

  signature_temp_rank <- exp_rank * gSet_pos_matrix
  signature_temp_rank <- signature_temp_rank[-which(rowSums(signature_temp_rank) == 0), ]
  signature_rank <- sortMatrixByCol(signature_temp_rank, decreasing=TRUE)

  score_for_oneset <- sapply(seq_len(ncol(exp_rank)), function(i){
    t4 <- indicatorInsideGeneSet_matrix[,i]
    t4[which(t4 != 0)] <- signature_rank[,i]
    stepCDFinGeneSet <- cumsum((abs(t4))^alpha)/ sum((abs(signature_rank[,i]))^alpha)
    t3 <- indicatorInsideGeneSet_matrix[,i]
    stepCDFoutGeneSet <- cumsum(!t3)/sum(!t3)
    walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
    return(sum(walkStat))
  })
  return(score_for_oneset)
}

