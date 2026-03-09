#' Generate a single simulated count matrix from a given graph (ZINB source model)
#'
#' Simulates one count matrix \eqn{X} associated with an input graph \code{net}.
#' The generator creates a matrix of latent sources \eqn{Y} with \eqn{m = p + p(p-1)/2}
#' columns, sampled from a zero-inflated negative binomial distribution, and then maps
#' these sources to \eqn{p} observed variables through a deterministic construction based
#' on the graph adjacency structure. An additional zero-inflated noise matrix \eqn{E} is added.
#'
#' The current implementation expects \code{net} to be an undirected graph with \code{p} vertices.
#' Internally, it uses the lower-triangular part of the adjacency matrix to encode edge structure.
#'
#' @param net An \code{igraph} object representing the ground-truth network. Must have \code{p} vertices.
#' @param seed Integer seed used to set the RNG state before simulating the matrix.
#' @param n Integer number of samples (rows) to generate.
#' @param p Integer number of variables / nodes (columns) in the generated matrix (graph size).
#' @param theta Positive numeric dispersion parameter passed to \code{emdbook::rzinbinom} as \code{size}.
#' @param means Numeric mean parameter passed to \code{emdbook::rzinbinom} as \code{mu} for the latent sources.
#'        In this implementation a single value is used for all latent-source draws.
#' @param zerop Numeric in \eqn{[0,1]} giving the zero-inflation probability \code{zprob} used in \code{emdbook::rzinbinom}.
#'
#' @return A numeric matrix \code{X} with \code{n} rows and \code{p} columns containing simulated counts.
#'
#' @examples
#' library(igraph)
#' set.seed(1)
#' p <- 6
#' g <- sample_gnp(p, 0.3, directed = FALSE)
#' X <- Generetwork_NB(
#'   net = g,
#'   seed = 123,
#'   n = 50,
#'   p = p,
#'   theta = 2,
#'   means = 1.5,
#'   zerop = 0.4
#' )
#' dim(X)
#'
#' @seealso \code{\link{k_net_datagen}} to generate multiple matrices from the same network.
#'
#' @export
#' @importFrom igraph as_adjacency_matrix
#' @importFrom emdbook rzinbinom
#' @importFrom utils combn
Generetwork_NB <- function(net, # ground truth
                           seed, # to have the same matrix if desired
                           n, # how many samples
                           p, # how many variables (graph size)
                           theta, # overdispersion of sources
                           means, # NB sources mean
                           zerop # zero percentage
                           )
  {
  # We'll have the graph realisation X as Y%*%B, so let's start from the dimension of matrix Y; an n row times p + p(p-1)/2 columns sources matrix
  bino <- choose(p, 2) # p(p-1)/2
  m <- p + bino # p + p(p-1)/2

  # We set the seed
  set.seed(seed)
  matseed <- .Random.seed

  # We can now generate Y matrix values:
  Y <- matrix(emdbook::rzinbinom(n = n * m,
                            mu = means, # means from park chosen 1.5
                            size = theta,
                            zprob = zerop),
              nrow = n,
              ncol = m)

  E <- matrix(emdbook::rzinbinom(n = n*p, mu = 0, size = 1e8, zprob = zerop), nrow = n, ncol = p)

  # We now got to build matrix B, this will be made by an identity matrix p times p, a permutation matrix and the vectorization of the upper triangular matrix of the graph adjacency matrix

  Ip <- diag(p)

  # To generate the permutation matrix let's start from all the possible two element combination
  combinations <- combn(1:p, 2)


  # And then we can get the permutation matrix
  permutation_matrix <- matrix(0, nrow = p, ncol = choose(p, 2))
  for (i in 1:ncol(permutation_matrix)) {
    permutation_matrix[combinations[1,i], i] <- 1
    permutation_matrix[combinations[2,i], i] <- 1
  }

  # Finally we retrieve the upper triangular vectorization of the adjacency matrix by row order

  # We exploit the fact that under diagonal elements of a symmetric matrix taken by column correspond to the upper one take by row
  adj_matrix <- igraph::as_adjacency_matrix(net, type = "lower", sparse = FALSE)

  # With these elements we can now build the fundamental vector giving our network structure to the data

  upper_tri_vector <- adj_matrix[lower.tri(adj_matrix)]
  upper_tri_matrix <- t(replicate(p, upper_tri_vector))

  # We cast the graph structure into the data by taking the product between the permutation matrix and the upper triangular matrix

  K <-permutation_matrix * upper_tri_matrix

  # Finally, to get matrix b, we vertically stuck K under the identity matrix [Ip]=pxp

  B <- t(cbind(Ip, K))

  # Our X matrix, simulating a realisation of the graph, possibly representing a scRNA-seq cluster is:

  X <- Y%*%B + E

  return(X = X)
}
#' Generate multiple simulated count matrices from the same network
#'
#' Convenience wrapper around \code{\link{Generetwork_NB}} to generate a list of \code{k}
#' matrices from the same ground-truth network. Each matrix uses a different seed (\code{seed + i})
#' and can use a different latent-source mean via \code{means[i]}.
#'
#' @param net An \code{igraph} object representing the ground-truth network. Must have \code{p} vertices.
#' @param k Integer number of matrices to generate.
#' @param seed Integer base seed. The i-th matrix is generated with seed \code{seed + i}.
#' @param n Integer number of samples (rows) in each generated matrix.
#' @param p Integer number of variables (columns) in each generated matrix.
#' @param theta Positive numeric dispersion parameter passed to \code{emdbook::rzinbinom} as \code{size}.
#' @param means Numeric vector of length \code{k}. The i-th element \code{means[i]} is passed to
#'        \code{\link{Generetwork_NB}} for the i-th matrix.
#' @param zerop Numeric in \eqn{[0,1]} giving the zero-inflation probability.
#'
#' @return A list of length \code{k}. Each element is a numeric matrix of dimension \code{n} x \code{p}.
#'
#' @examples
#' library(igraph)
#' p <- 8
#' g <- sample_gnp(p, 0.2, directed = FALSE)
#' Xlist <- k_net_datagen(
#'   net = g,
#'   k = 3,
#'   seed = 100,
#'   n = 40,
#'   p = p,
#'   theta = 2,
#'   means = c(1.2, 1.5, 1.8),
#'   zerop = 0.5
#' )
#' length(Xlist)
#' vapply(Xlist, dim, integer(2))
#'
#' @export
#' @seealso \code{\link{Generetwork_NB}}
k_net_datagen <- function(net,
                          k,
                          seed,
                          n,
                          p,
                          theta,
                          means,
                          zerop)
  {
  grpnetworks <- list()
  for (i in 1:k){
    grpnetworks[[i]] <- Generetwork_NB(net = net, seed = seed + i,
                                       n, p, theta, means = means[i], zerop = zerop)
    # net indica quale rete usare per la struttura dei dati
    # tramite seed, fissata una rete, la matrice risultante cambia
    # n è la sample size della matrice generata, ogni vertice produce n osservazioni
    # p è il numero di colonne ed equivale al numero di vertici del grafo
    # Theta controlla la dispersione mediante la formula var = mean + theta*mean^2
    # means controlla le medie delle sorgenti, vedi View(Generetwork_NB_grpnet)
  }
  return(grpnetworks)
}
#' Compute dropout rate (fraction of zeros) in a matrix
#'
#' Computes the proportion of entries equal to zero in a matrix \code{X}.
#'
#' @param X A numeric matrix (typically counts).
#'
#' @return A single numeric value in \eqn{[0,1]} equal to \code{sum(X == 0) / length(X)}.
#'
#' @examples
#' X <- matrix(c(0, 1, 0, 2), nrow = 2)
#' calculate_dropout(X)
#'
#' @export
calculate_dropout <- function(X) {
  zeroes <- sum(X == 0) / length(X)
  return(zeroes)
}
