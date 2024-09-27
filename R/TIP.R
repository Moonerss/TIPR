#' Run TIP analysis
#'
#' @description
#' Run TIP web server analysis
#'
#' @param expression expression matrix with gene in row and sample in column.
#' You'd better use \code{log2(tpm+1)} normalized data, or preprocess your data with \code{preprocess}
#' @param perm permutation times, default 100
#' @param seed random state,default 1234
#' @param verbose Gives information about each calculation step. Default: TRUE.
#' @param ncores the parallel number of cores to run analysis, default 1.
#'
#' @importFrom stringr str_subset str_split str_replace_all
#' @importFrom purrr map_chr
#'
#' @return return a matrix with TIP activity score
#'

TIP <- function(expression, perm = 100, seed = 1234, verbose = TRUE, ncores = 1) {

  # check parallel
  ncores <- as.integer(ncores)
  if (ncores > 1 & !requireNamespace("pbapply", quietly = TRUE)) {
    cli::cli_abort('parallel need `pbapply` package, please install it')
  }
  # set parallel
  if (ncores > 1) {
    if (.Platform$OS.type == "windows") {
      cli::cli_alert_warning(cli::col_yellow('Parallel is not supported on Windows'))
    }
  }

  gene_set_list <- TIPR::TIP_signature_symbol
  # intersect number of each signature sets
  geneset_intersect <- sapply(gene_set_list, function(i){
    num <- length(intersect(i, rownames(expression)))
    return(num)
  })

  # for step with positive and negative gene set,
  # if one direction of a step matched 0 genes,
  # the other direction will be also discarded
  if (verbose) cli::cli_alert_info('Filtering gene set ...')
  positive_geneset <- str_subset(names(gene_set_list), 'positive')
  negative_geneset <- str_subset(names(gene_set_list), 'negative')

  zero_intersect <- names(which(geneset_intersect == 0))
  positive_zero <- intersect(zero_intersect, positive_geneset)
  negative_zero <- intersect(zero_intersect, negative_geneset)

  if(length(positive_zero > 0) | length(negative_zero > 0)){
    full_zero_geneset <- c(positive_zero, negative_zero)
    delete_step <- full_zero_geneset %>% str_split('[.]') %>% purrr::map_chr(1) %>% unique()
    delete_pos_neg_signature <- str_subset(names(gene_set_list), pattern = paste0(delete_step, collapse = '|'))
    delete_signature <- c(delete_pos_neg_signature, zero_intersect) %>% unique()
  }else{
    delete_signature <- zero_intersect %>% unique()
  }

  if (length(delete_signature) > 0) {
    if (verbose) {
      delete_fun <- delete_signature %>%
        str_replace_all('.positive', '') %>%
        str_replace_all('.negative', '') %>%
        unique()
      cli::cli_alert_warning('Next {.val {length(delete_fun)}} signature{?s} will be filtered: {.val {delete_fun}}')
    }
    gene_set_list <- gene_set_list[-delete_signature]
  }

  # build random matrix, random gene
  if (verbose) cli::cli_alert_info('Random {.val {perm}} time{?s} ...')
  sample_number <- ncol(expression)
  if (ncores > 1) {
    perm_exp <- pbapply::pblapply(seq_len(perm), function(x) {
      multi_sample_permutation(expression, seed = seed)
    }, cl = ncores)
  } else {
    perm_exp <- lapply(seq_len(perm), function(x) {
      multi_sample_permutation(expression, seed = seed)
    })
  }

  perm_exp <- Reduce(cbind, perm_exp)
  perm_exp <- cbind(expression, perm_exp) # add raw expression

  rownames(perm_exp) <- rownames(expression)

  # get the superimposed matrix
  superposition_rank <- superpositionRank(matrix = perm_exp)
  # get rank by column of matrix 'perm_exp'
  perm_exp_rank <- rankMatrixByCol(superposition_rank = superposition_rank)

  # activity score of signature sets one by one
  if (verbose) cli::cli_alert_info('Calculate activity score ...')
  if (ncores > 1) {
    permutation_score <- t(pbapply::pbsapply(gene_set_list, function(sig){
      gSetIdx <- which(rownames(perm_exp) %in% sig)
      # vector of '0'and '1' show the position of signature genes in 'perm_exp'
      t1 <- numeric(nrow(perm_exp))
      t1[gSetIdx] <- 1
      geneSet_pos <- matrix(rep(t1, times = ncol(perm_exp)), ncol = ncol(perm_exp))
      # activity score of one signature set
      score_for_oneset <- ssGSEA_ES(
        exp_rank = perm_exp_rank,
        superposition_rank = superposition_rank,
        gSet_pos_matrix = geneSet_pos)

      return(score_for_oneset)
    }, cl = ncores))
  } else {
    permutation_score <- t(sapply(gene_set_list, function(sig){
      gSetIdx <- which(rownames(perm_exp) %in% sig)
      # vector of '0'and '1' show the position of signature genes in 'perm_exp'
      t1 <- numeric(nrow(perm_exp))
      t1[gSetIdx] <- 1
      geneSet_pos <- matrix(rep(t1, times = ncol(perm_exp)), ncol = ncol(perm_exp))
      # activity score of one signature set
      score_for_oneset <- ssGSEA_ES(
        exp_rank = perm_exp_rank,
        superposition_rank = superposition_rank,
        gSet_pos_matrix = geneSet_pos)

      return(score_for_oneset)
    }))
  }


  # normalize activity score
  if (verbose) cli::cli_alert_info('Normalize activity score ...')
  ## 1. whether activity score of real sample and permutation samples(N+1，2N+1 ... perm.times*N+1) are of the same sign,
  ## if different, set the score of permutation samples to 'NA'
  sample_num <- ncol(expression)
  real_sample_score <- permutation_score[, seq_len(sample_num)]
  for (i in seq_len(perm)) {
    one_perm_score <- permutation_score[, i*sample_num + seq_len(sample_num)]
    r <- which(real_sample_score * one_perm_score < 0)
    one_perm_score[r] <- NA
    permutation_score[, i*sample_num + seq_len(sample_num)] <- one_perm_score
  }

  ## 2. zscore normalization
  zscore_perm <- sapply(seq_len(sample_num), function(s){
    #activity score of one sample and it's permutation sample
    one_sample_perm_mat <- permutation_score[,((0:perm) * sample_num + s)]
    #z-score normalization
    zscore_one_sample <- apply(one_sample_perm_mat, 1, function(x) {
      if(length(which(!is.na(x))) < 2){
        zscore <- x[1]
      }else{
        zscore <- scale(x)[1]
      }
      return(zscore)
    })
    return(zscore_one_sample)
  })

  ## 3. positive and negative gene set
  pos_fun <- str_subset(rownames(zscore_perm), "positive")
  pos_idx <- match(pos_fun, rownames(zscore_perm))
  neg_fun <- str_replace_all(pos_fun, 'positive', 'negative')
  neg_idx <- match(neg_fun, rownames(zscore_perm))

  if(length(pos_fun) > 0){
    #change the sign of score of negative genes
    zscore_perm[neg_fun, ] <- zscore_perm[neg_fun, ]*(-1)
    direction <- seq_len(nrow(zscore_perm))
    direction[pos_idx] <- direction[neg_idx]
  }else{
    direction <- seq_len(nrow(zscore_perm))
  }
  #rowsum to conbine positive genes score and negative genes score
  step_names <- rownames(zscore_perm) %>%
    str_replace_all('.negative', '') %>%
    str_replace_all('.positive', '')
  normalized_score <- rowsum(zscore_perm, step_names)

  colnames(normalized_score) <- colnames(expression)
  if (verbose) cli::cli_alert_success('Done')
  return(normalized_score)
}
