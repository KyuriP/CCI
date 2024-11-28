#' FCI's skeleton discovery procedure adapted for CCI
#' @param suffStat List containing suffStat$data
#' @param indepTest The independence test
#' @param alpha Type I error rate for indepTest
#' @return List containing the skeleton data

IP_discovery <- function(suffStat, indepTest, alpha, p, max.cs = Inf, 
                         fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, ...) {
  
  # Start timing the process
  time_start <- proc.time()
  
  # # Match the JCI assumption
  # jci <- match.arg(jci)
  # 
  # # Handle JCI assumptions by setting fixed edges between context and system variables
  # if (jci == "123" && !is.null(contextVars)) {
  #   # Ensure fixedEdges is a square matrix
  #   if (is.null(fixedEdges)) {
  #     fixedEdges <- matrix(FALSE, p, p)
  #   }
  #   # Set fixed edges between context variables for JCI 123 assumption
  #   fixedEdges[contextVars, contextVars] <- TRUE
  # }
  # 
  # if (jci %in% c("1", "123") && !is.null(contextVars)) {
  #   # Direct edges between context and system variables (JCI 1 and 2)
  #   for (i in contextVars) {
  #     for (j in setdiff(seq_len(p), contextVars)) {
  #       fixedEdges[i, j] <- TRUE
  #     }
  #   }
  # }
  
  # Run the skeleton discovery using the modified skeleton function
  skel <- skeleton_new_jci(suffStat, indepTest, alpha, p = p, 
                           m.max = max.cs, 
                           fixedGaps = NULL, fixedEdges = NULL, verbose = verbose, ...)
  
  # Extract the graph and separation sets from the skeleton phase
  G_sk <- as(skel@graph, "matrix")
  sepset_sk <- skel@sepset
  
  # Measure time for the skeleton phase
  time_skel <- proc.time() - time_start
  
  # Perform the possible-d-separation (pdsep) step based on the skeleton
  pdsepRes <- pdsep(skel@graph, suffStat, indepTest, p = p, sepset = sepset_sk, 
                    alpha = alpha, pMax = skel@pMax, m.max = max.cs, 
                    fixedEdges = fixedEdges, verbose = verbose, ...)  # Pass fixedEdges here
  
  # Extract the final graph and separation sets after the pdsep phase
  G <- pdsepRes$G
  sepset <- pdsepRes$sepset
  
  # Measure time for the entire process
  time_pdsep <- proc.time() - time_start
  
  # Store the results in a list and return
  resIP <- list(
    G = G,
    G_sk = G_sk,
    sepset = sepset,
    sepset_sk = sepset_sk,
    time_skel = time_skel,
    time_pdsep = time_pdsep,
    skel = skel,
    pdsepRes = pdsepRes
  )
  
  return(resIP)
}

