#' FCI's skeleton discovery procedure adapted for CCI
#' @param suffStat List containing suffStat$data
#' @param indepTest The independence test
#' @param alpha Type I error rate for indepTest
#' @return List containing the skeleton data

IP_discovery <- function(suffStat, indepTest, alpha, p, max.cs = Inf, contextVars = NULL, jci = c("0", "1", "123"), verbose = FALSE) {
  
  # Start timing the process
  time_start <- proc.time()
  
  # Run the skeleton discovery using the modified skeleton function that incorporates JCI assumptions
  skel <- skeleton_new_jci(suffStat, indepTest, alpha, p = p, m.max = max.cs, contextVars = contextVars, jci = jci, verbose = verbose)
  
  # Extract the graph and separation sets from the skeleton phase
  G_sk <- as(skel@graph, "matrix")
  sepset_sk <- skel@sepset
  
  # Measure time for the skeleton phase
  time_skel <- proc.time() - time_start
  
  # Perform the possible-d-separation (pdsep) step based on the skeleton
  pdsepRes <- pdsep(skel@graph, suffStat, indepTest, p = p, sepset = sepset_sk, alpha = alpha, pMax = skel@pMax, m.max = max.cs, verbose = verbose)
  
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