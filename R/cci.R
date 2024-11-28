
#' Cyclic Causal Inference (CCI) algorithm.
#'
#' Discovers a partially oriented maximal almost ancestral graph (MAAG) of a directed graph G,
#' provided that the global Markov property and d-separation faithfulness holds according to G.
#'
#' @param suffStat list of sufficient statistics needed by the CI test. E.g., the data or covariance matrix
#' @param indepTest the CI test function
#' @param alpha the alpha value for the CI test
#' @param p the number of variables
#' @param skeleton_pre Pre-computed skeleton
#' @param rules A logical vector indicating which of the 7 orientation rules to fire
#' @return A list containing the partially oriented MAAG \code{paag} and statistic \code{Sta}
#' @export

library(kpcalg)

cci <- function(suffStat, indepTest, alpha, p, skeleton_pre = NULL, 
                rules = rep(TRUE, 7), verbose = FALSE, labels = NULL, 
                jci = c("0", "1", "12", "123"), contextVars = NULL, fixedGaps = NULL, fixedEdges = NULL) {
  
  # Match the JCI assumption
  jci <- match.arg(jci)
  
  # Validate and handle context variables
  if (!is.null(contextVars) && length(contextVars) > 0) {
    if (!is.numeric(contextVars) || any(contextVars != as.integer(contextVars)) || 
        min(contextVars) < 1 || max(contextVars) > p) {
      stop("contextVars must be a vector of integers between 1 and p.")
    }
  }
  
  # Initialize fixed edges based on JCI assumption 3
  fixedEdges <- NULL
  if (jci == "123" && length(contextVars) > 0) {
    fixedEdges <- matrix(FALSE, p, p)
    fixedEdges[contextVars, contextVars] <- TRUE
  }
  
  # Step 1: Discover skeleton using IP_discovery (already modified for JCI)
  if (verbose) {
    cat("Compute Skeleton with JCI\n========================\n")
  }
  
  pdsepRes <- IP_discovery(suffStat, indepTest = indepTest, alpha = alpha, p = p, 
                           fixedEdges = fixedEdges)
  
  G <- pdsepRes$G
  sepset <- pdsepRes$sepset
  pMax <- pdsepRes$pMax
  allPdsep <- pdsepRes$allPdsep
  n.edgetestsPD <- pdsepRes$n.edgetests
  max.ordPD <- pdsepRes$max.ord
  
  if (is.null(labels)) {
    labels <- as.character(seq_len(p))
  }
  
  tripleList <- NULL
  
  # Step 2: Orient v-structures (already considering unfounded triples)
  G <- v_struc(pag = G, sepset, unfVect = tripleList, verbose)
  
  # Step 3: Post V-Structure adjustment
  list_pre_v <- after_v_struc(pag = G, sepset, suffStat, indepTest, alpha, verbose = verbose)
  
  # Step 4: Additional D-separation logic
  sup_sepset <- extra_dsep(list_pre_v$G, suffStat, indepTest, sepset, alpha, verbose = verbose)
  
  # Step 5: Refine graph using JCI and sepset
  list_e <- step_5(list_pre_v$G, sepset, sup_sepset, suffStat, indepTest, alpha, verbose = verbose,
                   rules_used = list_pre_v$rules_used)
  
  # Step 6: Further update graph with separation sets
  list_f <- step_6(list_e$pag, sepset, sup_sepset, suffStat, indepTest,
                   alpha, verbose = verbose, list_e$rules_used)
  
  # Step 7: Final edge orientation with udag2pag4 (already modified for JCI)
  res <- udag2pag4(pag = list_f$pag, sepset, rules = rules, unfVect = tripleList, 
                   verbose = verbose,  contextVars = contextVars, jci = jci)
  
  # Ensure result is a numeric matrix and label it
  res$pag <- apply(res$pag, 2, as.numeric)
  dimnames(res$pag) <- list(labels, labels)
  
  # Return results with additional info
  return(list(maag = res$pag, pre_OR_res = list_f$pag, 
              rules_used = res$rules_used, sepset = sepset))
}
