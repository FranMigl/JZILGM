#' In this script we have all the accessory functions to run the method J-ZINB-LGM
#'
#' The "indexes" function extracts from the input list of matrices \code{Xs}, all with same number of \code{p}
#' variables, the necessary indexes to generate a sparse design matrix suitable for a grouped lasso regression
#'
#' @param Xs a list of count matrices, these are the datasets to integrate
#' @param p common genes of interest
#' @return a list of useful structural indexes for a design matrix
#' @keywords internal
#' @noRd
indexes <- function(Xs, p){
  stopifnot(is.list(Xs))
  if (!is.numeric(p) || length(p) != 1L || p <= 0) stop("p must be a positive integer.")
  p <- as.integer(p)

  k = length(Xs)

  # Basic validation: all matrices and all have p columns
  ncols <- vapply(Xs, ncol, integer(1))
  if (anyNA(ncols) || any(ncols != p)) {
    stop("All matrices in Xs must have exactly p columns.")
  }

  # Total number of rows and columns of the design matrix
  totrow <- sum(sapply(Xs, nrow))
  totcol <- sum(sapply(Xs, ncol))

  # Cumulative sums
  row_indices <- cumsum(sapply(Xs, nrow))
  col_indices <- cumsum(sapply(Xs, ncol))

  ##################################### Intercept matrix:

  pattern_matrix <- Matrix::Matrix(data = 0, nrow =  totrow, ncol = k, sparse = TRUE)
  for (j in 1:k) {
    pattern_matrix[row_indices[j] - nrow(Xs[[j]]) + 1:nrow(Xs[[j]]), j] <- 1
  }

  ##################################### Design Matrix without intercepts

  matricione <- Matrix::Matrix(data = 0, nrow = totrow, ncol = totcol, sparse = TRUE)
  # in "matricione" we're diagonally inserting the k matrices
  for (j in 1:k) {
    matricione[row_indices[j] - nrow(Xs[[j]]) + 1:nrow(Xs[[j]]),
               col_indices[j] - ncol(Xs[[j]]) + 1:ncol(Xs[[j]])] <- Xs[[j]]
  }


  # Final design matrix

  matricione <- cbind(pattern_matrix, matricione)


  nk <- sapply(Xs, nrow) # samples of each matrix
  ncsum = cumsum(nk) # cumulative number of samples
  inds <- Map(function(start, end) start:end, c(1, utils::head(ncsum, -1) + 1), ncsum)
  # k long list with to address row indexes per original matrix in the design matrix

  # grpnet requires group variables to be consecutives e.g. 111222333 etc.
    group = c(1:k, rep((k+1):(p+k), each = k))

  # Index that order the matrix according to the groups;

  ind <- as.numeric(c(1:k, sapply(1:(p), function(i){
    k + seq(from = i, to = k * p, by = p)
  })))

  # Same index withou intercepts
  ind_j <- as.numeric(sapply(1:(p), function(i){
    seq(from = i, to = k * p, by = p)
  }))

  # total number of column of the design matrix P

  P = ncol(matricione)

  output= list()
  output$matricione = matricione
  output$totrow = totrow
  output$totcol = totcol
  output$nk = nk
  output$ncsum = ncsum
  output$inds = inds
  output$group = group
  output$ind = ind
  output$ind_j = ind_j
  output$P = P
  return(output)
  # Useful indexes are:
  # totrow totcol: matrices indexes without intercepts
  # row_indices col_indices: dimensions increment for each matrix intercepts excluded
  # nk: non incremental matrices dimensions
  # ncsum: equal to row_indices (redundant)
  # inds: indexes to address each matrix in the design matrix
  # group: group indexes for grpnet
  # ind: to indicisize the matrix according to the group variable
  # ind_j: to fill coefficients matrix excluding intercepts
  # P: final variables dimension of the design matrix
  # n, p = number of rows (genes), and variables (genes)
}
#'
#' Compute lambda_max for grouped design blocks (internal)
#' @param X_list List of design matrices (each with same number of columns).
#' @return Numeric scalar lambda_max.
#' @keywords internal
#' @noRd
find_lammax_grp = function(X_list)
{
  # let's FAST compute T(X) %*% X
  Grammatrices <- lapply(X_list, crossprod)

  k = length(X_list)
  # Let's initialise the matrix containing the norms
  n = nrow(Grammatrices[[1]])

  norm_matrix <- matrix(0, nrow = n, ncol = n)

  # Let's compute the norms
  for (i in 1:n) {
    for (j in 1:n) {
      ijs <- sapply(Grammatrices, function(mat) mat[i, j])
      norm_matrix[i, j] <- sqrt(sum(ijs^2))
    }
  }
  # numsample <- norm(as.matrix(sapply(X_list, nrow)), type = "2")^2
  # numsample <- norm(as.matrix(sapply(X_list, nrow)), type = "2")
  numsample <- sqrt(k)*sum(sapply(X_list, nrow))*(length(X_list)^2) # best in simul regime
  # numsample <- mean(as.matrix(sapply(X_list, nrow)))
  rhomax = max(abs(norm_matrix[upper.tri(norm_matrix)]))
  #print(paste("RHOmax = ", rhomax, sep = ""))
  lammax = 1/numsample * rhomax
  #print(paste("Lammax: ", lammax, sep = ""))
  return(lammax)
}
#'
#' The next function combines grouped-lasso coefficients into joint norms (internal)
#'
#' Given the coefficient array from a grouped lasso regression, this function
#' computes a \code{p x p} matrix of joint coefficients by taking the
#' Euclidean (L2/Frobenius) norm across the \code{k} group-specific coefficients
#' for each variable pair, at a fixed penalization index.
#'
#' This is used internally to reconstruct adjacency matrices from
#' group-structured regression outputs.
#'
#' @param output_matrix A 3D numeric array of regression coefficients with
#'   dimensions \code{(k*p) x p x n_lambda}, where the first dimension stacks
#'   group-specific coefficients.
#' @param indgrp A list of length \code{p}; each element contains the row indices
#'   in \code{output_matrix} corresponding to one group of size \code{k}.
#' @param p Integer; number of variables (genes).
#' @param lambda Integer index selecting the penalization level along the third
#'   dimension of \code{output_matrix}.
#'
#' @return A \code{p x p} numeric matrix where entry \code{(i, j)} is the L2 norm
#'   of the \code{k} coefficients associated with variable pair \code{(i, j)}
#'   at penalization index \code{lam_idx}.
#'
#' @keywords internal
#' @noRd
per_lambda_normal <- function(output_matrix, indgrp, p, lambda) {
  # Initialize the result matrix for a single lambda
  norm_matrix_single <- matrix(0, nrow = p, ncol = p)

  # Loop over i and j to compute the norm for each element
  for (i in 1:p) {
    for (j in 1:p) {

      # Extract the corresponding coefficients from the output matrix for the current lambda
      coefficients <- output_matrix[indgrp[[i]],j,lambda]

      # Compute the norm of the extracted coefficients and store in the result matrix
      norm_matrix_single[i, j] <- sqrt(sum(coefficients^2))
    }
  }

  return(norm_matrix_single)
}
#' Negative binomial pmf (NB1 parametrization) (internal)
#'
#' Computes the negative binomial probability mass function under a mean/dispersion
#' parametrization using \code{mu} (mean) and \code{theta} (dispersion).
#' The implementation works on the log scale internally and exponentiates unless
#' \code{log = TRUE}.
#'
#' @param y Numeric vector of count values (assumed non-negative integers).
#' @param mu Numeric vector of negative binomial means (same length as \code{y} or recyclable).
#' @param theta Numeric dispersion parameter (same length as \code{y} or recyclable).
#' @param log Logical; if \code{TRUE}, return log-pmf values.
#'
#' @return A numeric vector of pmf (or log-pmf) values with the same length as \code{y}.
#'
#' @keywords internal
#' @noRd
dNBI = function(y, mu, theta, log = FALSE)
{
  density = lgamma(y + theta + 1e-10) - lgamma(y + 1) - lgamma(theta + 1e-10) +
    theta * (log(theta + 1e-10) - log(theta + mu + 1e-10)) +
    y * (log(mu + 1e-10) - log(theta + mu + 1e-10))

  if (log == FALSE) {density = exp(density)}

  return(density)
}
#' Negative binomial pmf (NB2 reparametrization) (internal)
#'
#' Wrapper around \code{dNBI} implementing the NB2 reparametrization where
#' \code{sigma = mu / theta}. Internally it calls \code{dNBI(y, mu, theta = mu/sigma)}.
#'
#' @param y Numeric vector of count values (assumed non-negative integers).
#' @param mu Numeric vector of means (same length as \code{y} or recyclable).
#' @param sigma Numeric overdispersion parameter (\code{sigma = mu/theta}).
#' @param log Logical; if \code{TRUE}, return log-pmf values.
#'
#' @return A numeric vector of pmf (or log-pmf) values with the same length as \code{y}.
#'
#' @keywords internal
#' @noRd
dNBII = function(y, mu, sigma, log = FALSE)
{
  density = dNBI(y, mu = mu, theta = mu / sigma, log = log)
  return(density)
}
#' Grouped ZINB/NB2 pseudo-objective (internal)
#'
#' Computes the pseudo-objective minimized during fitting under the NB2 parametrization.
#' The objective combines a (zero-inflated) negative binomial pseudo-loglikelihood term
#' and a grouped penalty computed from the grouped coefficient matrix \code{bvec}.
#'
#' This function expects the data to be partitioned across \code{k} datasets using
#' \code{inds} (row indices per dataset).
#'
#' @param y Numeric response vector of counts (concatenated across datasets).
#' @param weights Numeric vector of observation weights (same length as \code{y}).
#' @param prob Numeric vector of length \code{k}; probability of structural zeros per dataset.
#' @param bvec Numeric matrix of estimated grouped coefficients (as returned by the group-lasso solver).
#' @param mu Numeric vector of fitted means (same length as \code{y}).
#' @param sigma Numeric scalar or vector; NB2 overdispersion (passed to \code{dNBII}).
#' @param lambda Numeric scalar; current penalization level.
#' @param penalty.factor Numeric vector of group penalty weights (aligned to group indices).
#' @param posz Numeric/logical vector indicating zero positions in \code{y} (e.g. \code{y == 0}).
#' @param k Integer; number of datasets.
#' @param p Integer; number of variables/genes.
#' @param inds List of length \code{k}; each element contains row indices in \code{y} belonging to dataset \code{i}.
#'
#' @return Numeric scalar: objective value (pseudo-negative loglikelihood + penalty).
#'
#' @keywords internal
#' @noRd
grp_nb2_objective = function(y, weights, prob, bvec, mu, sigma = NULL, lambda, penalty.factor, posz, k, p, inds)
{
  if (is.null(weights)) stop("weights must be provided (non-NULL) here.")
  # Firstly we divide the beta results in groups, intercepts excluded
  # These will be the values to be l2 norm computed for the penalty
  bgrp <- lapply(seq(k + 1, dim(bvec)[1], by = k), function(start) {
    matrix(bvec[start:(start + k - 1), ])
  })

  # Iteratively computing penalty
  penalty = 0
  for (i in 1:p){
    penalty = penalty + penalty.factor[i + k]*norm(bgrp[[i]])
  }

  # actual likelihood computation
  pnl = 0
  for (i in 1:k) {
    pnl = pnl - sum(weights[inds[[i]]] * log(prob[i] * posz[inds[[i]]] + (1 - prob[i]) * dNBII(y = y[inds[[i]]], sigma = sigma, mu = mu[inds[[i]]],log = FALSE) + 1e-10))
  }
  return(pnl/k + lambda*penalty)
}
#' Overdispersion estimation under NB2 parametrization (internal)
#'
#' Estimates the NB2 overdispersion parameter \code{sigma} by maximizing a weighted
#' log-likelihood over \code{sigma} using \code{optimize()} on a fixed interval.
#' On warnings/errors, the function falls back to a moment-based expression.
#'
#' @param y Numeric response vector of counts.
#' @param mu Numeric vector of fitted means (same length as \code{y}).
#' @param weights Numeric vector of observation weights (same length as \code{y}).
#'
#' @return Numeric scalar: estimated overdispersion parameter \code{sigma}.
#' @importFrom stats optimize var
#'
#' @keywords internal
#' @noRd
sigma_ml = function(y, mu, weights = NULL)
{
  if (length(weights) != length(y)) stop("weights must have same length as y.")

  NB2_theta = function(sigma, mu, y, weights) {
    return(sum(weights * dNBII(y = y, sigma = sigma, mu = mu, log = TRUE)))
  }

  # start = c(0.01)
  sigma = tryCatch(
    {result = optimize(NB2_theta, y = y, mu = mu, weights = weights, interval = c(1e-6, 1000), maximum = TRUE)
    sigma = result$maximum
    return(sigma)
    },
    warning = function(w) {
      cat("Warning message:", conditionMessage(w), "\n")
      sigma = (mean(y)^2)/(var(y)-mean(y))
      print(paste("Sigma after warning :", sigma, sep = " "))
      return(sigma)
    }, error = function(e) {
      cat("Error message:", conditionMessage(e), "\n")
      sigma = (mean(y)^2)/(var(y)-mean(y))
      print(paste("Sigma after Error :", sigma, sep = " "))
      return(sigma)
    })

  # If an error or warning occurred, result will contain the value returned by the error or warning handler.
  # If no error or warning occurred, result will contain the result of the code block.

  # sigma = ifelse(sigma <= 5e-5, 0, sigma)
  return(sigma)
}
#' Construct a symmetric adjacency matrix from coefficients (internal)
#'
#' Converts a coefficient matrix into a binary, symmetric adjacency matrix using
#' an absolute-value threshold. Symmetrization is performed using either:
#' \itemize{
#'   \item \code{type = "AND"}: edge present if both \code{(i,j)} and \code{(j,i)} exceed threshold
#'   \item \code{type = "OR"}: edge present if either \code{(i,j)} or \code{(j,i)} exceeds threshold
#' }
#'
#' @param coef_mat Numeric matrix of coefficients.
#' @param thresh Numeric threshold for edge inclusion (default \code{1e-6}).
#' @param type Character; one of \code{"AND"} or \code{"OR"} controlling symmetrization.
#'
#' @return A numeric (0/1) adjacency matrix of the same dimension as \code{coef_mat}.
#' if OR then an edge is considered present when at least one of the symmetric values is above the threshold
#' @keywords internal
#' @noRd
hat_net = function(coef_mat, thresh = 1e-6, type = c("AND", "OR"))
{
  type = match.arg(type)

  tmp_mat = abs(coef_mat) > thresh

  if (type == "AND") {
    res_mat = (tmp_mat * t(tmp_mat)) * 1
  }

  if (type == "OR") {
    res_mat = (tmp_mat + t(tmp_mat) > 0) * 1
  }
  return(res_mat)
}
#' Compute edge-recovery metrics for an estimated network
#'
#' Compares an estimated network to a ground-truth network and returns confusion
#' counts and standard performance rates. The comparison is performed on the
#' upper triangle (excluding the diagonal), treating the networks as undirected
#' and counting each edge once.
#'
#' Both inputs are binarized internally: any nonzero entry is treated as an edge.
#'
#' @param net Ground-truth network, either an \code{igraph} object or a square
#'   adjacency matrix (numeric/logical). Nonzero entries are treated as edges.
#' @param hat Estimated network as a square adjacency matrix (numeric/logical).
#'   Nonzero entries are treated as edges.
#'
#' @return A named numeric vector with confusion counts (ZE,TN,FP,NZ,FN,TP),
#'   the number of edges in the estimate (\code{N.edges}), and rates
#'   (TPR,FNR,TNR,FPR,Precision,Balanced_Accuracy,F1).
#'
#' @examples
#' net <- matrix(0, 3, 3)
#' net[1,2] <- net[2,1] <- 1
#' net[2,3] <- net[3,2] <- 1
#' hat <- net
#' metrics(net, hat)
#'
#' @export
metrics <- function(net, hat) {
  # Accept igraph or adjacency matrix for net
  if (inherits(net, "igraph")) {
    net <- igraph::as_adjacency_matrix(net, type = "both", sparse = FALSE)
  } else {
    net <- as.matrix(net)
  }

  hat <- as.matrix(hat)

  if (!all(dim(net) == dim(hat))) stop("net and hat must have the same dimensions.")
  if (nrow(net) != ncol(net)) stop("net must be square.")

  # Binarize
  net_bin <- (net != 0)
  hat_bin <- (hat != 0)

  # Use only undirected unique edges: upper triangle, no diagonal
  ut <- upper.tri(net_bin, diag = FALSE)

  net_u <- net_bin[ut]
  hat_u <- hat_bin[ut]

  ZE <- sum(!net_u)           # true zeros in net
  NZ <- sum(net_u)            # true edges in net
  TN <- sum(!net_u & !hat_u)
  FP <- sum(!net_u &  hat_u)
  FN <- sum( net_u & !hat_u)
  TP <- sum( net_u &  hat_u)

  # Rates (guard division by zero)
  TPR <- if (NZ > 0) TP / NZ else NA_real_
  FNR <- if (NZ > 0) FN / NZ else NA_real_
  TNR <- if (ZE > 0) TN / ZE else NA_real_
  FPR <- if (ZE > 0) FP / ZE else NA_real_
  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Balanced_Accuracy <- if (is.na(TPR) || is.na(TNR)) NA_real_ else (TPR + TNR) / 2
  F1 <- if (is.na(Precision) || is.na(TPR) || (Precision + TPR) == 0) NA_real_ else 2 * (Precision * TPR) / (Precision + TPR)

  x <- c(ZE = ZE, TN = TN, FP = FP, NZ = NZ, FN = FN, TP = TP,
         `N.edges` = (TP + FP),
         TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR,
         Precision = Precision, Balanced_Accuracy = Balanced_Accuracy, F1 = F1)

  x
}
#' Helper for lambda bisection: fit network and count edges (internal)
#'
#' Runs \code{zinb_LGM_net_grp()} at a given \code{current_lambda} and returns
#' the fitted object along with the total number of edges across the estimated
#' networks (as returned in \code{tmp_net$hat_net}).
#'
#' This helper is used during bisection/selection of a lambda producing a target
#' sparsity level.
#'
#' @param current_lambda Numeric scalar; current lambda tested in the bisection process.
#' @param Xlist List of input count matrices.
#' @param mat_indexes Output of \code{indexes()}, containing the design matrix and indices.
#' @param nvar Integer; number of variables/genes.
#' @param sym Character; symmetrization rule passed to \code{zinb_LGM_net_grp()}.
#' @param theta Numeric; parameter forwarded to \code{zinb_LGM_net_grp()}.
#' @param thresh Numeric; threshold forwarded to \code{zinb_LGM_net_grp()}.
#' @param weights_mat Passed through to \code{zinb_LGM_net_grp()} (currently required to be \code{NULL}).
#' @param penalty_mat Passed through to \code{zinb_LGM_net_grp()} (currently required to be \code{NULL}).
#' @param nCores Integer; number of cores for parallel computation.
#' @param offset Optional offset passed to \code{zinb_LGM_net_grp()}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{edges}{Total number of edges across estimated networks.}
#'     \item{tmp_net}{The fitted object returned by \code{zinb_LGM_net_grp()}.}
#'   }
#'
#' @keywords internal
#' @noRd
adjust_lambda <- function(
    current_lambda,
    Xlist,
    mat_indexes,
    nvar,
    sym,
    theta,
    thresh,
    nCores,
    offset
) {

  tmp_net <- zinb_LGM_net_grp(
    X = mat_indexes$matricione,
    lambda = current_lambda,
    sym = sym,
    theta = theta,
    thresh = thresh,
    nCores = nCores,
    p = nvar,
    totrow = mat_indexes$totrow,
    totcol = mat_indexes$totcol,
    P = mat_indexes$P,
    k = length(Xlist),
    nk = mat_indexes$nk,
    inds = mat_indexes$inds,
    group = mat_indexes$group,
    ind = mat_indexes$ind,
    ind_j = mat_indexes$ind_j,
    offset = offset,
    verbose = 0
  )

  edges <- sum(unlist(lapply(tmp_net$hat_net, function(x) {
    igraph::gsize(igraph::graph_from_adjacency_matrix(methods::as(x, "CsparseMatrix"), mode = "undirected"))
  })))

  list(edges = edges, tmp_net = tmp_net)
}
