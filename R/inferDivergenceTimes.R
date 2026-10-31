###############################################################################
# Package dependencies for each model:
# - Library mgcv for inferring divergence times with GAM (Generalized Additive Model).
# - Library scam for inferring divergence times with SCAM (Shape Constrained Additive Model).
# - Library minpack.lm for inferring divergence times with exponential and power-law functions.
###############################################################################

##################################################
# Tree related functions.
##################################################

#' Leaves grouped by direct descendants of a given node
#'
#' Get all leaves descending from a given node, grouped by its direct descendants.
#'
#' @param node_id An integer corresponding to the node ID.
#' @param tree A [`ConvenientTblTree`] object.
#'
#' @return A named list mapping each direct descendant node ID to the IDs of its descendant leaves (a list of lists).
#'
#' @author Priscila Biller
getLeaves <- function(node_id, tree){
  children <- tree$node[(tree$parent == node_id) & (tree$parent != tree$node)] |> purrr::set_names()
  lapply(children, function(child_id) {tidytree::offspring(tree, child_id, tiponly=TRUE, self_include=TRUE)$node})
}

#' Cumulative branch lengths
#'
#' This method computes the cumulative sum of branch lengths along the path 
#' from the root to each node in the tree.
#'
#' @param tree A [`ConvenientTblTree`] object.
#' @param rootId The node ID of the root (column `node` in `tree`).
#'
#' @return The given tree with an additional `cumBranchLen` column containing the cumulative sum of branch lengths for each node.
#'
#' @author Priscila Biller
computeCumulativeBranchLengths <- function(tree, rootId=NULL){
  if(is.null(rootId)){
    rootId <- tree[is.na(tree$branch.length),]$node
  }
  tree[, "cumBranchLen"] <- 0.0
  stack <- list(list(node=rootId,pathSoFar=0.0))
  while (length(stack) > 0) {
    # Pop element.
    node  <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    # Update node's cumulative sum of branch lengths.
    tree$cumBranchLen[tree$node == node$node] <- node$pathSoFar
    # Add node's children to the stack.
    children <- tree$node[(tree$parent == node$node) & (tree$parent != tree$node)]
    for(child in children){
      stack[[length(stack) + 1]] <- list(node=child, pathSoFar=(node$pathSoFar+tree$branch.length[tree$node == child]))
    }
  }
  tree
}

##################################################
# TimeTree related functions.
# Functions used to get times from TimeTree.
##################################################

#' Gets all pairs of nodes with TimeTree info for a single node.
#'
#' Computes all pairs of leaves that could provide information about the divergence time of one ancestral node.
#'
#' @param children A named list of all (informative) leaves descending from a given node, grouped by its direct descendants.
#'                 The names of the elements are the node IDs of the direct descendants.
#'
#' @return A matrix in which each row corresponds to a pair of leaves.
#'
#' @importFrom utils combn
#'
#' @author Priscila Biller
getTimeTreePairs_onenode <- function(children){
  # Computes all possible pairwise combinations of direct descendants.
  # Usually an internal/ancestral node has only two descendants, but here the code 
  # is a bit more generic in case there is an unresolved node with more than two descendants.
  children <- children[lengths(children) > 0]
  if (length(children) < 2) {
    return(matrix(character(0),nrow = 0,ncol = 2,dimnames = list(NULL, c("leaf_1", "leaf_2"))))
  }
  child_pairs    <- combn(names(children), 2, simplify = FALSE)
  # For each pair of direct descendants, generates all pairs of leaves such that one leaf comes from each descendant.
  all_leaf_pairs <- lapply(child_pairs, function(child_pair) {
    leaves_right <- children[[child_pair[1]]]
    leaves_left  <- children[[child_pair[2]]]
    # Returns a matrix with all combinations for a given pair.
    expand.grid(leaves_right,leaves_left)
    })
  # Combine matrices from pairs of direct descendants into a single matrix.
  all_leaf_pairs <- do.call(rbind, all_leaf_pairs)
  colnames(all_leaf_pairs) <- c("leaf_1", "leaf_2")
  # Discard uninformative pairs (both leaves have the same NCBI id).
  all_leaf_pairs <- all_leaf_pairs[all_leaf_pairs[, "leaf_1"] != all_leaf_pairs[, "leaf_2"],,drop = FALSE]
  all_leaf_pairs
}

#' Gets all pairs of nodes with TimeTree info for all nodes.
#'
#' Compute all pairs of leaves that could provide information about the divergence time of an ancestral node.
#'
#' @param nodes A list in which each element corresponds to an ancestral node. 
#'              Each element was computed with the function `getLeaves` and consists of a named list 
#'              that maps the ID of each direct descendant node to the IDs of its descendant leaves.
#' @param tree A [`ConvenientTblTree`] object.
#' @param ncbi_ids_timetree A list of NCBI ids that appear in the TimeTree database.
#'
#' @return A named list that maps each ancestral node to pairs of leaves, identified by their node IDs. 
#'         These pairs of leaves are informative for estimating the divergence time of that ancestor.
#'
#' @author Priscila Biller
getTimeTreePairs <- function(nodes, tree, ncbi_ids_timetree){

  # For each ancestral node, filter its descendant leaves to retain only those present in the TimeTree database.
  childrenInTimeTree <- lapply(nodes, function(children_groups) {
    lapply(children_groups, function(children) { # `children` is a list of node ids (!= ncbi ids).
      children_ncbi_ids <- unique(tree$NCBI_id[match(children, tree$node)])
      children_ncbi_ids[!is.na(children_ncbi_ids) & (children_ncbi_ids %in% ncbi_ids_timetree)]
      })})

  # For each ancestral node, it computes all pairs of leaves that could provide information about its divergence time.
  lapply(childrenInTimeTree, getTimeTreePairs_onenode)
}

#' Gets the most recent common ancestor (MRCA) of two nodes.
#'
#' Computes the most recent common ancestor of two nodes, given their timelines from TimeTree.
#'
#' @param taxid1 An integer: The NCBI id of one taxon
#' @param taxid2 An integer: The NCBI id of the other taxon
#' @param timelines A named list mapping each node to its timeline from TimeTree. 
#'                  Each timeline contains all ancestors of the node up to the root of the Tree of Life.
#'                  A timeline is represented as a dataframe with one row per ancestor.
#'
#' @return The row in the dataframe corresponding to the most recent common ancestor of the two given nodes.
#'
#' @author Priscila Biller
getMRCA <- function(taxid1, taxid2, timelines){
  # Gets the TimeTree timelines of both nodes.
  timeline_1 <- timelines[[as.character(taxid1)]] |> dplyr::filter(divergenceTime > 0)
  timeline_2 <- timelines[[as.character(taxid2)]] |> dplyr::filter(divergenceTime > 0)
  # Computes the intersection of two timelines (= all common ancestors) and 
  # selects the row with the lowest divergence time, which corresponds to the most recent common ancestor (MRCA).
  dplyr::intersect(timeline_1,timeline_2) |> dplyr::slice_min(order_by=divergenceTime, n=1)
}

##################################################
# Different models to estimate divergence time.
# Functions used to extrapolate/interpolate times from TimeTree.
##################################################

#' Generalized Additive Model (GAM).
#' 
#' This function uses the relationship between cumulative branch length and 
#' time to infer the divergence times of ancestral nodes that have only the 
#' cumulative branch length.
#'
#' The relationship between cumulative branch length and time is inferred 
#' with GAM (Generalized Additive Model).
#' 
#' @param known Nodes with both divergence time (column `divergenceTime`) and cumulative branch length (column `cumBranchLen`).
#' @param unknown Nodes with only cumulative branch length (column `cumBranchLen`; column `divergenceTime` is `NA`).
#' @param x_grid Optional parameter with extra points beyond those with unknown divergence times. These can be used, for example, to plot or check accuracy over a wider range.
#'
#' @return A list with 2 types of predictions: 
#'         1. `pred_unknown`: predictions for nodes in the path with unknown divergence time;
#'         2. `pred_grid`: predictions for the extra points passed in `x_grid`.
#'
#' @importFrom mgcv gam s
#'
#' @author Priscila Biller
model_gam <- function(known, unknown, x_grid=NULL){
  # Fit data with GAM (Generalized Additive Model).
  fit_gam     <- gam(divergenceTime ~ s(cumBranchLen,k=4), data=known, method="REML")
  # Make predictions.
  pred_gam_unknown <- pmax(predict(fit_gam,newdata = data.frame(cumBranchLen = unknown$cumBranchLen)),0) # avoid negative predictions.
  pred_gam_grid    <- if(!is.null(x_grid)) pmax(predict(fit_gam,newdata = data.frame(cumBranchLen = x_grid)),0) else NULL
  list(pred_unknown=pred_gam_unknown, pred_grid=pred_gam_grid)
}

#' Shape Constrained Additive Model (SCAM)
#' 
#' This function uses the relationship between cumulative branch length and 
#' time to infer the divergence times of ancestral nodes that have only the 
#' cumulative branch length.
#'
#' The relationship between cumulative branch length and time is inferred 
#' with SCAM (Shape Constrained Additive Model). The constrain added here
#' enforces monotonic decrease.
#' 
#' @param known Nodes with both divergence time (column `divergenceTime`) and cumulative branch length (column `cumBranchLen`).
#' @param unknown Nodes with only cumulative branch length (column `cumBranchLen`; column `divergenceTime` is `NA`).
#' @param x_grid Optional parameter with extra points beyond those with unknown divergence times. These can be used, for example, to plot or check accuracy over a wider range.
#'
#' @return A list with 2 types of predictions: 
#'         1. `pred_unknown`: predictions for nodes in the path with unknown divergence time;
#'         2. `pred_grid`: predictions for the extra points passed in `x_grid`.
#'
#' @note This function is similar to `model_gamdec` (same principle but different results).
#'
#' @importFrom scam scam
#' @importFrom mgcv s
#' @importFrom stats Gamma
#'
#' @author Priscila Biller
model_scam <- function(known, unknown, x_grid=NULL){
  # Data for Gamma function needs to be strictly positive.
  known_positive <- subset(known, divergenceTime > 0)
  k <- 4
  # Fit data with SCAM.
  fit_scam <- scam(divergenceTime ~ s(cumBranchLen, bs="mpd", k=k),data = known_positive,family=Gamma(link="log"))
  # Make predictions.
  pred_scam_unknown <- pmax(expm1(predict(fit_scam, newdata = data.frame(cumBranchLen = unknown$cumBranchLen))),0) # avoid negative predictions.
  pred_scam_grid    <- if(!is.null(x_grid)) pmax(expm1(predict(fit_scam, newdata = data.frame(cumBranchLen = x_grid))),0) else NULL
  list(pred_unknown=pred_scam_unknown, pred_grid=pred_scam_grid)
}

#' Generalized Additive Model (GAM) with a decreasing monotonicity constraint.
#' 
#' This function uses the relationship between cumulative branch length and 
#' time to infer the divergence times of ancestral nodes that have only the 
#' cumulative branch length.
#'
#' The relationship between cumulative branch length and time is initially  
#' inferred with GAM (Generalized Additive Model), and later the constrain 
#' monotonic decrease is added.
#' 
#' @param known Nodes with both divergence time (column `divergenceTime`) and cumulative branch length (column `cumBranchLen`).
#' @param unknown Nodes with only cumulative branch length (column `cumBranchLen`; column `divergenceTime` is `NA`).
#' @param x_grid Optional parameter with extra points beyond those with unknown divergence times. These can be used, for example, to plot or check accuracy over a wider range.
#'
#' @return A list with 2 types of predictions: 
#'         1. `pred_unknown`: predictions for nodes in the path with unknown divergence time;
#'         2. `pred_grid`: predictions for the extra points passed in `x_grid`.
#'
#' @note This function is similar to `model_scam` (same principle but different results).
#'
#' @importFrom mgcv gam s smoothCon mono.con pcls Predict.matrix
#'
#' @author Priscila Biller
model_gamdec <- function(known, unknown, x_grid=NULL){

  # Set up the size of the basis functions/number of knots
  k <- 4
  # Fits the unconstrained GAM model, but gets the smoothness parameters that will be used later.
  unc <- gam(divergenceTime ~ s(cumBranchLen, k = k, bs = "cr"), data=known)

  # Creates the cubic spline basis functions of `x`
  # It returns an object containing the penalty matrix for the spline
  # among other things; see ?smooth.construct for description of each
  # element in the returned object
  sm <- smoothCon(s(cumBranchLen, k = k, bs = "cr"), known, knots = NULL)[[1]]

  # Gets the constraint matrix and constraint vector that imposes
  # linear constraints to enforce montonicity on a cubic regression spline
  # the key thing you need to change is `up`.
  # `up = TRUE` == increasing function
  # `up = FALSE` == decreasing function
  # `xp` is a vector of knot locations that we get back from smoothCon
  F <- mono.con(sm$xp, up = FALSE)   # get constraints: up = FALSE == Decreasing constraint!

  # Fill in G, the object pcsl needs to fit; this is just what `pcls` says it needs:
  # X is the model matrix (of the basis functions)
  # C is the identifiability constraints - no constraints needed here
  #   for the single smooth
  # sp are the smoothness parameters from the unconstrained GAM
  # p/xp are the knot locations again, but negated for a decreasing function
  # y is the response data
  # w are weights and this is fancy code for a vector of 1s of length(y)
  G <- list(X = sm$X, C = matrix(0,0,0), sp = unc$sp,
            p = -sm$xp, # note the - here! This is for decreasing fits!
            y = known$divergenceTime,
            w = known$divergenceTime*0+1)
  G$Ain <- F$A    # the monotonicity constraint matrix
  G$bin <- F$b    # the monotonicity constraint vector, both from mono.con
  G$S <- sm$S     # the penalty matrix for the cubic spline
  G$off <- 0      # location of offsets in the penalty matrix

  # Fit the constrained spline.
  fit_coefficients <- pcls(G)  # fit spline (using s.p. from unconstrained fit)

  # Make predictions.
  pred_gamdec_unknown <- pmax(as.vector(Predict.matrix(sm, data.frame(cumBranchLen = unknown$cumBranchLen)) %*% fit_coefficients),0)
  pred_gamdec_grid    <- if(!is.null(x_grid)) pmax(as.vector(Predict.matrix(sm, data.frame(cumBranchLen = x_grid)) %*% fit_coefficients),0) else NULL
  list(pred_unknown=pred_gamdec_unknown, pred_grid=pred_gamdec_grid)
}

#' Isotonic regression.
#' 
#' This function uses the relationship between cumulative branch length and 
#' time to infer the divergence times of ancestral nodes that have only the 
#' cumulative branch length.
#'
#' The relationship between cumulative branch length and time is inferred
#' with Isotonic regression. Isotonic regression is the technique of fitting 
#' a free-form line to a sequence of observations such that the fitted line 
#' is non-decreasing (or non-increasing) everywhere, and lies as close to 
#' the observations as possible.
#' 
#' @param known Nodes with both divergence time (column `divergenceTime`) and cumulative branch length (column `cumBranchLen`).
#' @param unknown Nodes with only cumulative branch length (column `cumBranchLen`; column `divergenceTime` is `NA`).
#' @param x_grid Optional parameter with extra points beyond those with unknown divergence times. These can be used, for example, to plot or check accuracy over a wider range.
#'
#' @return A list with 2 types of predictions: 
#'         1. `pred_unknown`: predictions for nodes in the path with unknown divergence time;
#'         2. `pred_grid`: predictions for the extra points passed in `x_grid`.
#'
#' @importFrom stats isoreg approx
#'
#' @author Priscila Biller
model_iso <- function(known, unknown, x_grid=NULL){
  # Sort observations by x coords (i.e. cumulative branch lenght).
  known <- known[order(known$cumBranchLen), ]
  # Fit isotonic regression.
  fit_iso <- isoreg(x=known$cumBranchLen, y=-known$divergenceTime)
  # Prediction function for the isotonic curve
  predict_iso <- function(x_new) {
    pred <- approx(x=known$cumBranchLen, y=fit_iso$yf, xout=x_new, method="constant", f=0, rule=2)$y
    # Convert back from -divergenceTime to divergenceTime
    as.vector(-pred)
  }
  # Make predictions.
  pred_iso_unknown <- predict_iso(unknown$cumBranchLen)
  pred_iso_grid    <- if(!is.null(x_grid)) predict_iso(x_grid) else NULL
  list(pred_unknown=pred_iso_unknown, pred_grid=pred_iso_grid)
}

#' Linear interpolation.
#' 
#' This function uses the relationship between cumulative branch length and 
#' time to infer the divergence times of ancestral nodes that have only the 
#' cumulative branch length.
#'
#' The relationship between cumulative branch length and time is inferred
#' with Linear interpolation. 
#' 
#' @param known Nodes with both divergence time (column `divergenceTime`) and cumulative branch length (column `cumBranchLen`).
#' @param unknown Nodes with only cumulative branch length (column `cumBranchLen`; column `divergenceTime` is `NA`).
#' @param x_grid Optional parameter with extra points beyond those with unknown divergence times. These can be used, for example, to plot or check accuracy over a wider range.
#'
#' @return A list with 2 types of predictions: 
#'         1. `pred_unknown`: predictions for nodes in the path with unknown divergence time;
#'         2. `pred_grid`: predictions for the extra points passed in `x_grid`.
#'
#' @importFrom stats approx
#'
#' @author Priscila Biller
model_lin <- function(known, unknown, x_grid=NULL){
  # Sort observations by x coords (i.e. cumulative branch lenght).
  known <- known[order(known$cumBranchLen), ]
  # Make predictions.
  pred_lin_unknown <- approx(x=known$cumBranchLen, y=known$divergenceTime, xout=unknown$cumBranchLen, method="linear", rule=2)$y
  pred_lin_grid    <- if(!is.null(x_grid)) approx(x=known$cumBranchLen, y=known$divergenceTime, xout=x_grid, method="linear", rule=2)$y else NULL
  list(pred_unknown=pred_lin_unknown, pred_grid=pred_lin_grid)
}

#' Exponential model.
#' 
#' This function uses the relationship between cumulative branch length and 
#' time to infer the divergence times of ancestral nodes that have only the 
#' cumulative branch length.
#'
#' The relationship between cumulative branch length and time is assumed to 
#' follow a decreasing exponential curve.
#' 
#' @param known Nodes with both divergence time (column `divergenceTime`) and cumulative branch length (column `cumBranchLen`).
#' @param unknown Nodes with only cumulative branch length (column `cumBranchLen`; column `divergenceTime` is `NA`).
#' @param x_grid Optional parameter with extra points beyond those with unknown divergence times. These can be used, for example, to plot or check accuracy over a wider range.
#'
#' @return A list with 2 types of predictions: 
#'         1. `pred_unknown`: predictions for nodes in the path with unknown divergence time;
#'         2. `pred_grid`: predictions for the extra points passed in `x_grid`.
#'
#' @importFrom minpack.lm nlsLM
#' @importFrom stats median
#'
#' @author Priscila Biller
model_exp <- function(known, unknown, x_grid=NULL){
  # Fit exponential model.
  fit_exp <- nlsLM(divergenceTime ~ c + a * exp(-b * cumBranchLen), data=known,
    start = list(
      c = min(known$divergenceTime, na.rm=TRUE),
      a = max(known$divergenceTime, na.rm=TRUE) -
          min(known$divergenceTime, na.rm=TRUE),
      b = 1 / median(known$cumBranchLen, na.rm=TRUE)
    ),
    lower = c(c = 0, a = 0, b = 0),
    upper = c(c = Inf, a = Inf, b = Inf)
  )
  # Make predictions.
  pred_exp_unknown <- as.vector(predict(fit_exp, newdata=data.frame(cumBranchLen=unknown$cumBranchLen)))
  pred_exp_grid    <- if(!is.null(x_grid)) as.vector(predict(fit_exp, newdata=data.frame(cumBranchLen=x_grid))) else NULL
  list(pred_unknown=pred_exp_unknown, pred_grid=pred_exp_grid)
}

#' Power-law model.
#' 
#' This function uses the relationship between cumulative branch length and 
#' time to infer the divergence times of ancestral nodes that have only the 
#' cumulative branch length.
#'
#' The relationship between cumulative branch length and time is assumed to 
#' follow a decreasing power-law curve.
#' 
#' @param known Nodes with both divergence time (column `divergenceTime`) and cumulative branch length (column `cumBranchLen`).
#' @param unknown Nodes with only cumulative branch length (column `cumBranchLen`; column `divergenceTime` is `NA`).
#' @param x_grid Optional parameter with extra points beyond those with unknown divergence times. These can be used, for example, to plot or check accuracy over a wider range.
#'
#' @return A list with 2 types of predictions: 
#'         1. `pred_unknown`: predictions for nodes in the path with unknown divergence time;
#'         2. `pred_grid`: predictions for the extra points passed in `x_grid`.
#'
#' @importFrom minpack.lm nlsLM
#' @importFrom stats deviance
#'
#' @author Priscila Biller
model_power <- function(known, unknown, x_grid=NULL){
  # Try different start points.
  starts <- expand.grid(d=c(0.05, 0.1, 0.2), p=c(0.5, 1, 2))

  # Fit power-law model.
  # Parameters from the model:
  # c: curve value as cumBranchLen gets large;
  # a: curve height
  # d: behavior near cumBranchLen = 0;
  # p: curve shape
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    d0 <- starts$d[i]
    p0 <- starts$p[i]
    c0 <- min(known$divergenceTime)
    x0 <- known$cumBranchLen[which.min(known$cumBranchLen)]
    y0 <- known$divergenceTime[which.min(known$cumBranchLen)]
    a0 <- max(1e-6, (y0 - c0) * (x0 + d0)^p0)
    try(nlsLM(
      divergenceTime ~ c + a / (cumBranchLen + d)^p,
      data = known,
      start = list(c = c0, a = a0, d = d0, p = p0),
      lower = c(c = 0, a = 0, d = 1e-6, p = 0.05),
      upper = c(c = Inf, a = Inf, d = Inf, p = 5)
    ), silent = TRUE)
  })

  # Take the fits that did not fail.
  ok <- vapply(fits, inherits, logical(1), "nls")
  if (!any(ok)) {return(NA)}

  # Select the fit with the lowest residual.
  rss <- vapply(fits[ok], deviance, numeric(1))
  fit_power <- fits[ok][[which.min(rss)]]

  # Make predictions.
  pred_power_unknown <- as.vector(predict(fit_power, newdata=data.frame(cumBranchLen=unknown$cumBranchLen)))
  pred_power_grid    <- if(!is.null(x_grid)) as.vector(predict(fit_power, newdata=data.frame(cumBranchLen=x_grid))) else NULL
  list(pred_unknown=pred_power_unknown, pred_grid=pred_power_grid)
}

##################################################
# Plotting functions.
##################################################

#' Plot relation between divergence time and cumulative branch length for a clade.
#' 
#' This function saves a plot that shows the relationship between 
#' divergence time (taken from TimeTree) and cumulative branch length
#' for the whole clade.
#'
#' @param tree A [`ConvenientTblTree`] object.
#' @param cladeName The name of the analyzed clade.
#' @param fileFormat File format for saving the plot. Either `svg` or `png` options are available.
#'
#' @importFrom grDevices svg png dev.off
#'
#' @author Priscila Biller
plotDivergenceTimesFromTimeTree <- function(tree, cladeName, fileFormat="svg"){
  # Plot known divergence times from Time tree.
  if (fileFormat == "svg"){
    svg(filename=paste0(cladeName, "_obs.svg"), width=1800, height=1200)  
  } else {
    png(filename=paste0(cladeName, "_obs.png"), width=1800, height=1200, res=200) 
  }
  ancestors <- tree$node[tree$isTip == FALSE] |> purrr::set_names()
  observed  <- tree[(tree$node %in% ancestors) & (tree$isTimeExtrapolated == FALSE),]
  unknown   <- tree[(tree$node %in% ancestors) & (tree$isTimeExtrapolated == TRUE),]
  plot(observed$cumBranchLen, observed$divergenceTime, pch=19, col="black", xlab="Cumulative branch length", 
        ylab="Divergence time",main = paste0(cladeName, " - Timetree data: known=",nrow(observed),", unknown=",nrow(unknown)))
  dev.off()
}

#' Plot stats about each path from tip to root.
#' 
#' This function saves a plot that shows, for each path from tip to root:
#' 1. how many ancestral nodes have information from TimeTree;
#' 2. how many ancestral nodes do not have information from TimeTree because all tips belong to the same species;
#' 3. how many ancestral nodes do not have information from TimeTree, but tips are different species.
#' 
#' @param tree A [`ConvenientTblTree`] object.
#' @param pathEstimates Named list of dataframes, where each dataframe corresponds to a path from tip to root and its results.
#' @param cladeName The name of the analyzed clade.
#' @param fileFormat File format for saving the plot. Either `svg` or `png` options are available.
#'
#' @importFrom grDevices svg png dev.off
#' @importFrom graphics barplot text
#'
#' @author Priscila Biller
plotDivergenceTimesPerPath <- function(tree, pathEstimates, cladeName, fileFormat="svg"){
  # Join all dataframes in a single one.
  pathEstimates_all <- dplyr::bind_rows(pathEstimates,.id = "node_ref")
  # Create a new column to be used in the bar plot.
  pathEstimates_all <- pathEstimates_all |> dplyr::mutate(bartype = dplyr::case_when(
                                        (isTimeExtrapolated == FALSE) ~ "Not extrapolated",
                                        ((isTimeExtrapolated == TRUE) & (isSingleSpecies == FALSE)) ~ "Extrapolated, not single species",
                                        ((isTimeExtrapolated == TRUE) & (isSingleSpecies == TRUE)) ~ "Extrapolated, single species"
                                      ))
  bartype_levels <- c("Not extrapolated", "Extrapolated, not single species", "Extrapolated, single species")

  # Count how many known (isTimeExtrapolated == FALSE) and unknown  (isTimeExtrapolated == TRUE) 
  # divergence times there are for each path. Ancestral nodes of single species are put in a separated category.
  ancestors <- tree$node[tree$isTip == FALSE] |> purrr::set_names()
  counts    <- with(pathEstimates_all[pathEstimates_all$node %in% ancestors,],table(node_ref, factor(bartype, levels=bartype_levels)))
  
  # Sort paths by their proportion of information from TimeTree (isTimeExtrapolated == FALSE); paths with more information come first.
  counts <- counts[, bartype_levels, drop=FALSE]
  timeTreeInfoPerPath <- counts[, "Not extrapolated"] / rowSums(counts)
  nbNodesPerPath      <- rowSums(counts)

  sort_order <- order(timeTreeInfoPerPath, nbNodesPerPath, decreasing = TRUE, na.last = TRUE)
  counts_sorted <- counts[sort_order, ]

  if (fileFormat == "svg"){
    svg(filename=paste0(cladeName, "_obs_perpath.svg"), width=1800, height=1200)
  } else {
    png(filename=paste0(cladeName, "_obs_perpath.png"), width=1800, height=1200, res=200)
  }

  bar_positions <- barplot(t(counts_sorted),beside = FALSE,col = c("steelblue", "goldenrod", "tomato"),xlab = "Path",ylab = "Number of ancestral nodes",legend.text = bartype_levels, args.legend = list(x = "topright",bty = "n"))

  # Add labels inside each stacked bar.
  cumulative <- apply(counts_sorted, 1, cumsum)
  for (node_i in seq_len(nrow(counts_sorted))) {
    segment_bottom <- 0.0
    for (bar_j in seq_along(bartype_levels)) {
      value <- counts_sorted[node_i, bar_j]
      if (value > 0) {
        text(x=bar_positions[node_i], y=segment_bottom+value/2, labels=value, col="white", cex=0.8)
      }
      segment_bottom <- segment_bottom + value
    }
  }
  dev.off()
}

#' Plot the relation between distance from tips and lack of information from TimeTree.
#' 
#' @param unknownTypes List where each element corresponds to a path from tip to root. 
#'                     Each element in the list contains the distances from tips of all nodes 
#'                     without TimeTree information for the given path.
#'                     A distance equals to `1` means that the internal node is the direct parent of a leaf, 
#'                     whereas a higher distance `n` means that the node is closer to the root.
#' @param cladeName The name of the analyzed clade.
#' @param fileFormat File format for saving the plot. Either `svg` or `png` options are available.
#'
#' @importFrom grDevices svg png dev.off
#' @importFrom graphics points lines axis
#'
#' @author Priscila Biller
plotUnknownDistFromTips <- function(unknownTypes, cladeName, fileFormat="svg"){
  x_values <- rep(seq_along(unknownTypes),lengths(unknownTypes))
  y_values <- unlist(unknownTypes,use.names = FALSE)

  if (fileFormat == "svg"){
    svg(filename=paste0(cladeName, "_unk_tipdist.svg"), width=1800, height=1200)
  } else {
    png(filename=paste0(cladeName, "_unk_tipdist.png"), width=1800, height=1200, res=200)
  }

  # Set up an empty plot
  plot(x_values, y_values, type="n", xaxt="n", xlim=c(0.5, length(unknownTypes) + 0.5), ylim =range(y_values), xlab="Internal node with unknown time", ylab="Distance from leaves")

  # Add points and vertical lines for each list element
  for (i in seq_along(unknownTypes)) {
    y <- unknownTypes[[i]]
    x <- rep(i, length(y))
    points(x,y,pch = 19)
    # Connect points with the same x-coordinate
    lines(x,y,lwd = 2)
  }

  # Add labels for list elements
  axis(side=1, at=seq_along(unknownTypes), labels=names(unknownTypes))
  dev.off()
}

##################################################
# Functions to infer divergence time.
##################################################

#' Divergence time based on TimeTree 
#' 
#' Infer the divergence time of one node based on the TimeTree information.
#'
#' This method estimates the divergence time of a node given "informative" 
#' pairs of descendant leaves.
#'
#' An informative pair of leaves consists of two leaves, each belonging to a different 
#' direct descendant of the node, for which TimeTree data are available.
#'
#' For each pair, the method identifies the most recent common ancestor (MRCA) and retrieves 
#' its estimated divergence time. Ideally, all pairs should have the same MRCA (i.e., the given node) 
#' and yield the same divergence time estimate. If they do not, this indicates that the topology 
#' of the phylogeny from TimeTree differs from that inferred by the ScrambledTreeBuilder package.
#'
#' @param node_pairs A matrix with two columns (`leaf_1` and `leaf_2`) computed with the function `getTimeTreePairs`.
#'                   Each row in the matrix corresponds to a pair of informative leaves, identified by their node ID.
#' @param timelines A named list mapping each node to its timeline from TimeTree. 
#'                  Each timeline contains all ancestors of the node up to the root of the Tree of Life.
#'                  A timeline is represented as a dataframe with one row per ancestor.
#'
#' @return A data frame with one row per informative pair of leaves. 
#'         Each row is derived from the `timelines` object and contains the node IDs of the paired leaves, 
#'         along with information about their most recent common ancestor (MRCA), including its divergence time.
#'
#' @author Priscila Biller
inferDivTimesFromTimeTree <- function(node_pairs, timelines){
  # Gets the MRCA for all "informative" pairs of leaves descending from one ancestor. 
  # An informative pair consists of one leaf from one direct descendant and another leaf from a different direct descendant.
  # It is expected that all pairs point to the same MRCA (= ancestral node being evaluated). 
  # If one or more pairs have different MRCA, there is a mismatch between the TimeTree phylogeny and the tree inferred by the ScrambledTreeBuilder package.
  mrcas <- lapply(seq_len(nrow(node_pairs)), function(i) {
    leaf_1 <- node_pairs$leaf_1[i]
    leaf_2 <- node_pairs$leaf_2[i]
    mrca <- getMRCA(leaf_1,leaf_2,timelines = timelines)
    # Add two columns to the dataframe to identify which pair of leaves gives a certain divergence time estimate.
    mrca$leaf_1 <- leaf_1
    mrca$leaf_2 <- leaf_2
    mrca
    })
  # Joins all dataframes from the leaf pairs, creating one dataframe per ancestral node.
  mrcas <- do.call(rbind, mrcas)
  mrcas
}

#' Divergence time based on extrapolation
#'
#' This method infers the divergence time for nodes lacking them based on the divergence times from other nodes in the tree.
#'
#' The method first finds a function describing the relationship between cumulative branch length and divergence time.
#' It then uses this function to estimate divergence times for nodes that lack them, given their cumulative branch lengths.
#'
#' As the relationship between cumulative branch length and divergence time is not necessarily linear, the function uses 
#' a Generalized Additive Model (GAM) to infer this relationship.
#'
#' @param pathToRoot A [`ConvenientTblTree`] object containing all nodes along the path from a tip to the root of the tree.
#' @param cladeName An optional parameter that describes the dataset.
#' @param savePlots Boolean flag. If `TRUE`, the function saves `.png` and `.svg` files containing 
#'                  a plot with cumulative branch length on the x-axis and divergence time on the y-axis. 
#'                  The plot shows the fitted prediction function, predicted values, and observed data.
#' @param verbose Boolean flag (optional; default value is `FALSE`). 
#'                If `TRUE`, the function prints messages about the computation and summary statistics for the dataset.
#' @param expectedVals An optional dataframe with columns `cumBranchLen` and `divergenceTime` used for plotting.
#' @param models  An optional parameter defining which models will be used for extrapolating the data from TimeTree.
#'                Notice that only models with sufficient data are included in the results; the amount of information required varies by model. 
#'                Available options: `t_gam` (Generalized Additive Model),`t_scam` (Shape Constrained Additive Model),
#'                `t_gamdec` (Generalized Additive Model with a decreasing monotonicity constraint),
#'                `t_iso` (Isotonic regression), `t_lin` (Linear interpolation),
#'                `t_exp` (Exponential model),`t_power` (Power-law model).
#'
#' @return A named list with the following:
#'         1. `models` :  a vector of the models actually used to extrapolate the data from TimeTree. 
#'                        A model passed via the input parameter `models` is omitted if there is 
#'                        insufficient data to use it. `NULL` is returned if no model was used.
#'         2. `path` :  The [`ConvenientTblTree`] object received as input with divergence-time estimates assigned 
#'                      to ancestral nodes that lacked them.
#'         3. `known` : Number of points in the path with information (i.e., ancestral nodes with both cumulative branch length and divergence time).
#'         4. `unknown` : Number of points in the path without information (i.e., ancestral nodes with only cumulative branch length).
#'         5. `known_unique_cumblen` : Number of distinct cumulative branch lengths.
#'         6. `known_unique_divtime` : Number of distinct divergence times.
#'
#' @note  See functions `model_gam` (option `t_gam`), `model_scam` (option `t_scam`), `model_gamdec` (`t_gamdec`), 
#'        `model_iso` (option `t_iso`), `model_lin` (option `t_lin`), `model_exp` (option `t_exp`), 
#'        `model_power` (option `t_power`) for more details on each model.
#'
#' @importFrom grDevices png dev.off
#' @importFrom graphics lines points legend
#'
#' @author Priscila Biller
inferDivTimesFromExtrapolation <- function(pathToRoot, cladeName=NULL, savePlots=FALSE, verbose=FALSE,
      expectedVals=NULL, models=c("t_gam","t_scam","t_gamdec","t_iso","t_lin","t_exp","t_power")){

  # Sort nodes by cumulative branch length.
  pathToRoot     <- pathToRoot[order(pathToRoot$cumBranchLen), ]

  known          <- pathToRoot[!is.na(pathToRoot$divergenceTime), ]
  known          <- known[order(known$cumBranchLen), ]
  idx_unknown    <- is.na(pathToRoot$divergenceTime) & !is.na(pathToRoot$cumBranchLen) & (pathToRoot$cumBranchLen > 0)
  unknown        <- pathToRoot[idx_unknown, ]

  if (verbose) print(paste0("[",cladeName,"] Divergence times in a path from tip to root: known=",nrow(known),"(unique: cumblen=", length(unique(known$cumBranchLen)), "; divtime=", length(unique(known$divergenceTime)), ") / unknown=", nrow(unknown)))

  # Check if there are at least 2 points with known divergence time, and 1 point without:
  if ((length(unique(known$cumBranchLen)) < 2) || (length(unique(known$divergenceTime)) < 2) || (nrow(unknown) < 1)) {
    return(list(models=NULL, path=pathToRoot, 
                known=nrow(known), unknown=nrow(unknown),
                known_unique_cumblen=length(unique(known$cumBranchLen)), known_unique_divtime=length(unique(known$divergenceTime))))
  }

  # Different models to estimate divergence time.
  prediction_models_info <- list(
    t_gam    = list(func=model_gam,    color="dodgerblue",      minNbPoints=5),
    t_scam   = list(func=model_scam,   color="darkolivegreen4", minNbPoints=5),
    t_gamdec = list(func=model_gamdec, color="darkorange2",     minNbPoints=5),
    t_iso    = list(func=model_iso,    color="gold3",           minNbPoints=2),
    t_lin    = list(func=model_lin,    color="yellow2",         minNbPoints=2),
    t_exp    = list(func=model_exp,    color="deeppink",        minNbPoints=4),
    t_power  = list(func=model_power,  color="purple",          minNbPoints=4)
  )

  # Make predictions with different models.
  models_used      <- character()
  prediction_grids <- list()

  # Grid for plotting predictions.
  x_grid <- seq(min(c(known$cumBranchLen, pathToRoot$cumBranchLen), na.rm = TRUE), max(c(known$cumBranchLen, pathToRoot$cumBranchLen), na.rm=TRUE), length.out=300)
  known_range   <- range(known$cumBranchLen, na.rm = TRUE)
  inside_known  <- x_grid >= known_range[1] & x_grid <= known_range[2]
  outside_known <- x_grid < known_range[1] | x_grid > known_range[2]
  for (model_name in models) {
    model_info        <- prediction_models_info[[model_name]]
    model_function    <- model_info$func
    model_minNbPoints <- model_info$minNbPoints

    # Extrapolate divergence times only if at least `minPoints` distinct points with known divergence time exist. 
    # The number of minimum points depends on the method; for example, 4 is the minimum amount of information needed for `gam` (see parameter `k` in `gam`).
    if ((length(unique(known$cumBranchLen)) < model_minNbPoints) || (length(unique(known$divergenceTime)) < model_minNbPoints)) {
      next
    }

    # Estimate divergence times.
    if (!model_name %in% names(pathToRoot)) {pathToRoot[[model_name]] <- NA}
    models_used       <- c(models_used, model_name)
    model_predictions <- model_function(known=known, unknown=unknown, x_grid=x_grid)

    # Save predictions.
    pathToRoot[[model_name]][idx_unknown] <- model_predictions$pred_unknown
    prediction_grids[[model_name]] <- list(pred_grid=model_predictions$pred_grid, pred_unknown=model_predictions$pred_unknown)
  }

  # Plot known observations and predictions.
  if ((savePlots) && (length(models_used) > 0)){
    plotLabel <- paste0("node",pathToRoot$node[!pathToRoot$node %in% pathToRoot$parent])
    plotLabel <- ifelse(is.null(cladeName),    plotLabel, paste0(cladeName,"_",plotLabel))   
    plotLabel <- ifelse(is.null(expectedVals), plotLabel, paste0(plotLabel,"_groundtruth"))   
     
    png(filename=paste0(plotLabel, "_extrapolation.png"), width=1800, height=1200, res=200)

    # Combine prediction grids from all models
    pred_grids_all <- unlist(lapply(prediction_grids, function(model_grid) {model_grid$pred_grid}), use.names=FALSE)
    y_values <- c(known$divergenceTime, if (!is.null(expectedVals)) expectedVals$divergenceTime, pred_grids_all)
    y_values <- y_values[is.finite(y_values)]

    x_limits <- range(pathToRoot$cumBranchLen, na.rm=TRUE)
    y_limits <- range(y_values, na.rm=TRUE)

    # Plot observed values.
    plot(known$cumBranchLen,known$divergenceTime, pch=19, col="black", cex=1.5, xlim=x_limits, ylim=y_limits,
      xlab="Cumulative branch length", ylab="Divergence time", 
      main=paste0(plotLabel,", known=",nrow(known),"[ unique: cumblen=", length(unique(known$cumBranchLen)), "; divtime=", length(unique(known$divergenceTime)), "] / unk=", nrow(unknown))
      )
    # Plot predictions.
    colors      <- character()
    for (model_name in models_used) {
      model_info        <- prediction_models_info[[model_name]]
      model_color       <- model_info$color
      model_predictions <- prediction_grids[[model_name]]
      colors            <- c(colors, model_color)
      lines(x_grid[outside_known], model_predictions$pred_grid[outside_known], col=model_color, lwd=2, lty=2)  # Dashed line for extrapolation.
      lines(x_grid[inside_known],  model_predictions$pred_grid[inside_known],  col=model_color, lwd=2)         # Solid line for interpolation.
      points(unknown$cumBranchLen, model_predictions$pred_unknown,             col=model_color, pch=19, cex=1) # Dots for predictions.
    }
    # Plot expected values.
    if(!is.null(expectedVals)){
      points(expectedVals$cumBranchLen,expectedVals$divergenceTime, pch=4, col="black", lwd=2, cex=1.5)  
    }
    # Add legend.
    legend("topright",legend=models_used, col=colors, lwd=2, lty=1, bty="n", title="Model")
    dev.off()
  }
  list(models=models_used, path=pathToRoot, 
    known=nrow(known), unknown=nrow(unknown),
    known_unique_cumblen=length(unique(known$cumBranchLen)), known_unique_divtime=length(unique(known$divergenceTime)))
}

#' Estimate divergence time for ancestral nodes
#'
#' This method infers the divergence time of ancestral nodes based on the information from TimeTree when available, and
#' by extrapolating this information to other nodes without this information when possible.
#'
#' @param tree A [`ConvenientTblTree`] object.
#' @param timelines An optional named list of data frames, where each name corresponds to an NCBI ID. 
#'                  Each element in the list is a data frame containing the timeline from the specified 
#'                  taxon to the last common ancestor of all cellular organisms. 
#'                  Each row of the dataframe represents one ancestral node.
#'                  This list is computed by this function and returned in `timeTreeTimelines`, 
#'                  but because computation can take time, a previously computed list can also 
#'                  be given as an optional argument.
#' @param cladeName An optional parameter with the name of the analyzed clade.
#' @param savePlots Boolean flag (optional; default value is `FALSE`). 
#'                  If `TRUE`, the function saves `.png` and `.svg` files containing 
#'                  several plots (e.g. relationship between cumulative branch length and divergence time,
#'                  stats about the nodes for each path from tip to root, etc.) 
#' @param verbose Boolean flag (optional; default value is `FALSE`). 
#'                If `TRUE`, the function prints messages about the computation and summary statistics for the dataset.
#' @param models  An optional parameter defining which models will be used for extrapolating the data from TimeTree.
#'                Models must be listed in order of preference. The first model is preferred for setting the divergence time.
#'                If there is insufficient information to use a model, the next model is tried, and so on.
#'                Available options: `t_gam` (Generalized Additive Model),`t_scam` (Shape Constrained Additive Model),
#'                `t_gamdec` (Generalized Additive Model with a decreasing monotonicity constraint),
#'                `t_iso` (Isotonic regression), `t_lin` (Linear interpolation),
#'                `t_exp` (Exponential model),`t_power` (Power-law model).
#'
#' @return  A named list with the following information:
#'          1. `tree`:  A [`ConvenientTblTree`] object with the following additional columns: 
#'                      `divergenceTime` (the estimated divergence time for the node, either coming from TimeTree or extrapolated from TimeTree information), 
#'                      `isTimeExtrapolated` (`FALSE` if the divergence time comes directly from TimeTree, `TRUE` if it is extrapolated from TimeTree),
#'                      `cumBranchLen` (the sum of the edge lengths from the root to the node), 
#'                      `isSingleSpecies` (`TRUE` if all tips share the same `NCBI_id`; false otherwise),
#'                      plus additional columns, one per model used, named after the model and containing the divergence times inferred by that model.
#'          2. `timeTreeTimelines`: A list of data frames, where each name corresponds to an NCBI ID (from column `NCBI_id`). 
#'                                  Each element in the list is a data frame containing the timeline from the specified 
#'                                  taxon to the last common ancestor of all cellular organisms. 
#'                                  Each row of the dataframe represents one ancestral node.
#'          3. `timeTreeDivergenceTimes`: A named list of data frames, where each element corresponds to an ancestral node 
#'                                        in the tree for which a divergence time is available from TimeTree (names correspond to column `node`).
#'                                        Each element has a dataframe (computed with the function `inferDivTimesFromTimeTree`), 
#'                                        with one row per informative pair of leaves and the divergence time inferred from that pair.
#'                                        Notice that an ancestral node may have multiple pairs of leaves that yield different divergence times,
#'                                        indicating a possible disagreement between the TimeTree phylogeny and the phylogeny from ScrambledTreeBuilder.
#'          4. `pathsLeafToRoot`: All paths from tips to root used for the divergence time estimation, 
#'                                sorted by the amount of information from TimeTree they have.
#'
#' @note  See functions `model_gam` (option `t_gam`), `model_scam` (option `t_scam`), `model_gamdec` (`t_gamdec`), 
#'        `model_iso` (option `t_iso`), `model_lin` (option `t_lin`), `model_exp` (option `t_exp`), 
#'        `model_power` (option `t_power`) for more details on each model.
#'
#' @importFrom stats runif setNames
#'
#' @author Priscila Biller
inferDivergenceTimes <- function(tree, timelines=NULL, cladeName=NULL, savePlots=FALSE, verbose=FALSE, 
                                models=c("t_gamdec","t_lin","t_scam","t_exp","t_gam","t_iso")) {

  # Add new columns to the dataframe.
  tree[, "divergenceTime"]     <- ifelse(tree$isTip, 0, NA)
  tree[, "isTimeExtrapolated"] <- FALSE
  tree[, "cumBranchLen"]       <- 0.0
  root <- tree[is.na(tree$branch.length),]
  tree <- computeCumulativeBranchLengths(tree,rootId=root$node) 

  # Check which nodes are MRCA of a single species (i.e. all tips have the same NCBI id).
  tree[, "isSingleSpecies"] <- purrr::map_lgl(tree$node, function(node) {length(unique(tidytree::offspring(tree, node, tiponly=TRUE, self_include=TRUE)$NCBI_id)) == 1})

  # Add one column per model. Each column hold the estimate from that model.
  for (model in models) {
    tree[[model]] <- NA
  }

  # For each node with a NCBI id (should be leaves in most cases), try to fetch its timeline in the TimeTree database.
  if (is.null(timelines)) {
    if(verbose) print(paste("- Fetching data from TimeTree database..."))
    timelines <- tree$NCBI_id[!is.na(tree$NCBI_id)] |> unique() |> purrr::set_names() |> purrr::map(getTimeline)  
    timelines <- timelines[!sapply(timelines, function(timeline){length(timeline) == 1 && is.na(timeline)})]
  }
  
  cladeName <- ifelse(is.null(cladeName), runif(1, min=1, max=100), cladeName) 
  if(verbose) print(paste("- Clade ", cladeName, " has ", length(tree$node[tree$isTip == TRUE]), " leaves and ", length(timelines), " nodes with timelines."))
  
  # For each ancestor, it computes all pairs of leaves that could provide information about its divergence time.
  ancestors     <- tree$node[tree$isTip == FALSE] |> purrr::set_names()
  timeTreePairs <- ancestors |> lapply(getLeaves, tree=tree) |> getTimeTreePairs(tree=tree, ncbi_ids_timetree=names(timelines))

  # Counts how many pairs of informative leaves were found per ancestral node.
  timeTreePairsCounts <- setNames(vapply(timeTreePairs, nrow, integer(1)),names(timeTreePairs))
  timeTreePairsCounts <- sort(timeTreePairsCounts, decreasing = TRUE)

  # Step 1: go through all ancestors that have at least one child with TimeTree information on the left and one on the right.
  divergenceTimes <- timeTreePairs[names(timeTreePairsCounts[timeTreePairsCounts > 0])] |> lapply(inferDivTimesFromTimeTree, timelines=timelines)
  # Count the number of distinct divergence time estimates for each ancestral node.
  divergenceTimes_counts <- lapply(divergenceTimes, function(df) {table(df$divergenceTime)})

  nodes_TimeTree_singleTime <- 0
  nodes_TimeTree_manyTimes  <- 0
  nodes_TimeTree_singleTime_diffts <- length(unique(unlist(lapply(divergenceTimes_counts[lengths(divergenceTimes_counts) == 1], names), use.names = FALSE)))
  for (node in names(divergenceTimes_counts)) {
    # Check for tree conflicts between TimeTree and ScrambledTreeBuilder (i.e. if there is more than one distinct divergence time estimate for a node).
    counts_node <- table(divergenceTimes[[node]]$divergenceTime)
    row_idx     <- which(tree$node == as.numeric(node))
    if (length(counts_node) > 1) {
      # Print a warning for problematic nodes.
      print(paste("WARNING! Node ", node, " has multiple divergence time estimates. TimeTree info will be ignored and time will be extrapolated."))
      if(verbose) print(counts_node)
      nodes_TimeTree_manyTimes     <- nodes_TimeTree_manyTimes  + 1
      tree$divergenceTime[row_idx] <- NA
    } else {
      nodes_TimeTree_singleTime    <- nodes_TimeTree_singleTime + 1
      # Saves the maximum divergence time among the inferred values.
      # This value should match the value returned by the TimeTree API call named `MRCA`, 
      # accessible through: http://timetree.temple.edu/api/mrca/id/[ncbi ids separated by "+"]
      tree$divergenceTime[row_idx] <- max(divergenceTimes[[node]]$divergenceTime)
    }
  }

  # Step 2: Go through all ancestors that do not have informative pairs of leaves 
  # (i.e. no children with TimeTree information in either left or right subtree).
  tree[, "isTimeExtrapolated"] <- (is.na(tree$divergenceTime) & (tree$cumBranchLen > 0) & (tree$isTip == FALSE))
  leaves          <- tree$node[tree$isTip == TRUE]
  leafParents     <- tree$parent[match(leaves, tree$node)] # Goal is to remove duplicated paths (leaves with the same parent will have the same ancestors).
  pathsLeafToRoot <- leaves[!duplicated(leafParents)] |> purrr::set_names() |> lapply(function(leaf){dplyr::bind_rows(tree[tree$node == leaf, ], tidytree::ancestor(tree,leaf))})

  # Sort all paths from leaves to root by the amount of informative nodes they have.
  # Paths with more nodes containing TimeTree divergence times come first.
  pathsLeafToRoot <- pathsLeafToRoot[order(sapply(pathsLeafToRoot, function(path){length(unique(na.omit(path$divergenceTime)))}), decreasing=TRUE)]
  for (i in seq_along(pathsLeafToRoot)) {
    pathToRoot <- pathsLeafToRoot[[i]]
    # Check if path contains nodes with multiple times from TimeTree.
    expectedVals  <- NULL
    nodesTimeTree <- pathToRoot$node[is.na(pathToRoot$divergenceTime) & (as.character(pathToRoot$node) %in% names(divergenceTimes))]
    if (length(nodesTimeTree) > 0) {
      expectedVals <- data.frame(cumBranchLen = numeric(), divergenceTime = numeric())
      for (node in nodesTimeTree) {
        cumblen <- tree$cumBranchLen[tree$node == node]
        node_ts <- unique(divergenceTimes[[as.character(node)]]$divergenceTime)
        expectedVals <- rbind(expectedVals,data.frame(cumBranchLen = rep(cumblen, length(node_ts)),divergenceTime = node_ts))
      }
    }
    # Infer divergence times for nodes that are missing that information.
    extrapRes  <- inferDivTimesFromExtrapolation(pathToRoot=pathToRoot,cladeName=cladeName,savePlots=savePlots,expectedVals=expectedVals,models=models,verbose=verbose)
    if(!is.null(extrapRes$models)){
      pathToRoot <- extrapRes$path
      idx <- match(tree$node, pathToRoot$node)
      # Update tree with new estimates only if nodes in the tree do not have already an estimate.
      # This is done because paths are sorted by the amount of information from TimeTree they have;
      # thus paths processed first will have more reliable estimates, since the estimates are based
      # on more information available.
      update_idxs <- is.na(tree$divergenceTime) & !is.na(idx)
      for (model_name in extrapRes$models) {
        estimates    <- pathToRoot[[model_name]][idx]
        has_estimate <- update_idxs & !is.na(estimates)
        # Update tree with new estimates.
        tree[[model_name]][has_estimate] <- estimates[has_estimate]
        update_idxs_model <- is.na(tree$divergenceTime) & !is.na(estimates)
        tree$divergenceTime[update_idxs_model]     <- tree[[model_name]][update_idxs_model]
      }
      # Update path with new estimates.
      pathsLeafToRoot[[i]] <- pathToRoot
    }
  }

  # Saves plots and prints stats.
  if(savePlots){
    # Plot known divergence times from Time tree.
    plotDivergenceTimesFromTimeTree(tree, cladeName, fileFormat="svg")
    plotDivergenceTimesFromTimeTree(tree, cladeName, fileFormat="png")

    # Plot how many known and unknown divergence times there are per path from leaf to root.
    plotDivergenceTimesPerPath(tree, pathsLeafToRoot, cladeName, fileFormat="svg")
    plotDivergenceTimesPerPath(tree, pathsLeafToRoot, cladeName, fileFormat="png")
  }

  # Print stats on internal nodes.
  nodes_Extrapolation_singleSpecies <- sum((tree$isTimeExtrapolated == TRUE) & (tree$isSingleSpecies == TRUE)  & (tree$isTip == FALSE), na.rm=TRUE)
  nodes_Extrapolation_manySpecies   <- sum((tree$isTimeExtrapolated == TRUE) & (tree$isSingleSpecies == FALSE) & (tree$isTip == FALSE), na.rm=TRUE)
  nodes_Extrapolation_manySpecies   <- nodes_Extrapolation_manySpecies - nodes_TimeTree_manyTimes # Removing the nodes with TimeTree info that were also extrapolated.
  nodes_noEstimate_singleSpecies    <- sum(is.na(tree$divergenceTime) & (tree$isSingleSpecies == TRUE) & (tree$isTip == FALSE), na.rm=TRUE)
  nodes_noEstimate_manySpecies      <- sum(is.na(tree$divergenceTime) & (tree$isSingleSpecies == FALSE) & (tree$isTip == FALSE), na.rm=TRUE)
  if(verbose){
    print("Stats on internal nodes:")
    print(paste0("- Nodes w/ TimeTree and a single time    = ", nodes_TimeTree_singleTime, "; Number of distinct ts = ", nodes_TimeTree_singleTime_diffts))
    print(paste0("- Nodes w/ TimeTree and multiple times   = ", nodes_TimeTree_manyTimes))
    print(paste0("- Nodes no TimeTree and a single species = ", nodes_Extrapolation_singleSpecies))
    print(paste0("- Nodes no TimeTree and multiple species = ", nodes_Extrapolation_manySpecies))
    print(paste0("- Nodes no estimate and a single species = ", nodes_noEstimate_singleSpecies))
    print(paste0("- Nodes no estimate and multiple species = ", nodes_noEstimate_manySpecies))
  }

  list(tree=tree, timeTreeTimelines=timelines, timeTreeDivergenceTimes=divergenceTimes, pathsLeafToRoot=pathsLeafToRoot)
}


##################################################
# Functions to check divergence time inference.
##################################################

#' Performs checks on the estimated divergence times.
#' 
#' This function performs two different checks:
#' 1. Checks where ancestral nodes without TimeTree information occur most 
#'    often along the tree path: (a) near the root, (b) in the middle, or (c) near the tips.
#' 2. Checks the accuracy of different models for extrapolating the TimeTree information,
#'    by "erasing" some known points and comparing the estimates with the true values.
#'
#' @param tree A [`ConvenientTblTree`] object.
#' @param timelines An optional named list of data frames, where each name corresponds to an NCBI ID. 
#'                  Each element in the list is a data frame containing the timeline from the specified 
#'                  taxon to the last common ancestor of all cellular organisms. 
#'                  Each row of the dataframe represents one ancestral node.
#'                  This list is computed by this function and returned in `timeTreeTimelines`, 
#'                  but because computation can take time, a previously computed list can also 
#'                  be given as an optional argument.
#' @param cladeName An optional parameter with the name of the analyzed clade.
#' @param models  An optional parameter defining which models will be used for extrapolating the data from TimeTree.
#'                Models must be listed in order of preference. The first model is preferred for setting the divergence time.
#'                If there is insufficient information to use a model, the next model is tried, and so on.
#'                Available options: `t_gam` (Generalized Additive Model),`t_scam` (Shape Constrained Additive Model),
#'                `t_gamdec` (Generalized Additive Model with a decreasing monotonicity constraint),
#'                `t_iso` (Isotonic regression), `t_lin` (Linear interpolation),
#'                `t_exp` (Exponential model),`t_power` (Power-law model).
#' 
#' @return A list with two items: 
#'        1. `estimatesRes`: the output of the function `inferDivergenceTimes` (see its documentation for more details);
#'        2. `accuracySummary`: a dataframe with one row per tested model summarizing the accuracy of each tested model (e.g., mean, median, maximum, and minimum error).
#'
#' @importFrom stats median
#'
#' @author Priscila Biller
checkDivergenceTimeEstimates <- function(tree, timelines=NULL, cladeName=NULL, models=c("t_gamdec","t_lin","t_scam","t_exp","t_gam","t_iso","t_power")) {

  estimatesRes <- inferDivergenceTimes(tree, timelines=timelines, cladeName=cladeName, savePlots=TRUE, verbose=TRUE, models=models)

  tree            <- estimatesRes$tree
  timelines       <- estimatesRes$timelines
  pathsLeafToRoot <- estimatesRes$pathsLeafToRoot

  ###########################
  # Check 1: Check which type of unknown points appears more often:
  # (1) Unknown points near the root of the tree;
  # (2) Unknown points in the middle of the path;
  # (3) Unknown points in the end of the path.
  # This influences the choice of the model to extrapolate the known data.
  unknownDistFromTips <- lapply(pathsLeafToRoot, function(pathToRoot) {
                        known       <- pathToRoot[(!pathToRoot$isTimeExtrapolated) & (pathToRoot$isTip == FALSE), ]
                        # Check if path has information missing and at which point in the tree the information is missing.
                        if(nrow(known) < nrow(pathToRoot)){
                          # Sort nodes by cumulative branch length.
                          pathToRoot  <- pathToRoot[order(pathToRoot$cumBranchLen, decreasing=TRUE), ]
                          # 1 means parent of a leaf, n means root.
                          idx_unknown <- which(pathToRoot$isTimeExtrapolated & (pathToRoot$isTip == FALSE) & !is.na(pathToRoot$cumBranchLen) & (pathToRoot$cumBranchLen > 0))
                        } else {
                          NULL
                        }
                      })
  unknownDistFromTips <- Filter(function(x) !is.null(x),unknownDistFromTips)
  if(length(unknownDistFromTips) > 0){
    plotUnknownDistFromTips(unknownDistFromTips, cladeName, fileFormat="png")  
  }
  print("Location of unknown nodes (1=close to leaf; n=root)")
  print(unknownDistFromTips)

  ###########################
  # Check 2: Check accuracy of different models:
  # "Erase" some known points to check how accurate different methods are.
  accuracyRes <- data.frame()
  for(pathToRoot in pathsLeafToRoot){
    # Candidate nodes to be erased: (1) nodes with divergence time directly from TimeTree (not extrapolated), (2) not a tip, (3) not a root.
    known      <- pathToRoot[(!pathToRoot$isTimeExtrapolated) & (pathToRoot$isTip == FALSE) & (pathToRoot$cumBranchLen > 0), ]
    min_points_inference <- 4
    # Check if path has enough information:
    # At least 5 known points: 4 for the inference, and at least 1 point is going to be erased.
    if((length(unique(known$cumBranchLen)) > min_points_inference) & (length(unique(known$divergenceTime)) > min_points_inference)){ # (nrow(known) > (nrow(pathToRoot)*0.25)) & 
      # print(paste("Known=",nrow(known)," / Total=",nrow(pathToRoot)))
      # Sort nodes by cumulative branch length and remove the tip (does not need its divergence time estimated).
      pathToRoot <- pathToRoot[!pathToRoot$isTip, ]
      pathToRoot <- pathToRoot[order(pathToRoot$cumBranchLen), ]
      # Erases randomly 1, 2, ..., n rows (nodes can be at any point in the tree).
      max_erased <- min(length(unique(known$cumBranchLen)),length(unique(known$divergenceTime)))-min_points_inference  # At least 4 nodes w/information left for the inference.
      rep_cnt    <- 3 # Nb. repetitions for each nb. of erased points
      if(max_erased > 0){
        comb_samples <- expand.grid(erased=seq_len(max_erased), repetition=seq_len(rep_cnt))
        known_idx    <- which((!pathToRoot$isTimeExtrapolated) & (pathToRoot$isTip == FALSE) & (pathToRoot$cumBranchLen > 0))
        pathToRoot_samples <- lapply(seq_len(nrow(comb_samples)), function(row_idx) {
          pathToRoot[sample(known_idx, size=comb_samples$erased[row_idx], replace=FALSE),,drop = FALSE]
        })
        names(pathToRoot_samples) <- sprintf("n%d_m%d",comb_samples$erased,comb_samples$repetition)
        for (i in seq_along(pathToRoot_samples)) {
          # Holds all divergence times that will be erased (ground truth).
          testName <- names(pathToRoot_samples)[[i]]
          nodes_to_be_erased  <- pathToRoot_samples[[i]]
          # Erases divergence times.
          pathToRoot_erased <- pathToRoot
          pathToRoot_erased$divergenceTime[match(nodes_to_be_erased$node, pathToRoot_erased$node)] <- NA
          # Infer erased times.
          extrapRes <- inferDivTimesFromExtrapolation(pathToRoot=pathToRoot_erased,cladeName=paste0(cladeName,"_",testName),savePlots=TRUE,expectedVals=nodes_to_be_erased,models=models)
          # Compute accuracy.
          pathToRoot_res  <- extrapRes$path
          groundTruth_idx <- match(nodes_to_be_erased$node, pathToRoot_res$node)
          accuracyRes_i   <- data.frame(testName=testName, node=nodes_to_be_erased$node, true_divergenceTime=nodes_to_be_erased$divergenceTime)
          for(model in models){
            accuracyRes_i[[model]] <- abs(pathToRoot_res[[model]][groundTruth_idx]-accuracyRes_i$true_divergenceTime)
          }
          accuracyRes <- rbind(accuracyRes, accuracyRes_i)
        }
      }
    }
  }
  accuracySummary <- NA
  if(nrow(accuracyRes) > 0){
    accuracySummary <- data.frame(model = models, 
      sum = sapply(accuracyRes[models], sum, na.rm=TRUE),
      average = sapply(accuracyRes[models], mean, na.rm=TRUE),
      median = sapply(accuracyRes[models], median, na.rm=TRUE),
      q25 = sapply(accuracyRes[models], quantile, probs=0.25, na.rm=TRUE),
      q75 = sapply(accuracyRes[models], quantile, probs=0.75, na.rm=TRUE),
      max = sapply(accuracyRes[models], max, na.rm=TRUE),
      min = sapply(accuracyRes[models], min, na.rm=TRUE),
      row.names = NULL
    )
    accuracySummary <- accuracySummary[order(accuracySummary$median), ]
    print("Accuracy: Summary")
    print(accuracySummary) 
  }
  list(estimatesRes=estimatesRes,accuracySummary=accuracySummary)
}
