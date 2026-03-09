#' The first following function just run the method over temporary lambdas, if a grid of lambda is not provided.
#' These provisory lambdas then allow to obtain a more refined grid over which the method is then run
#' @param Xlist list of input data matrices, with the same variables
#' @param lambda optional grid of penalization values
#' @param conperc integer percentage of connectivity between 1 and 100 of the fully connected network of size \eqn{p(p-1)/2}. Default to 10.
#' @param sym how to symmetrize the output matrices, AND implies we consider an edge present if both the symmetrical values are above a certain threshold, OR consider an edge present if at least one of the symmetrical values is above the threshold
#' @param offset this is the model offset parameter, must be the total number of counts for a certain cell above all the data matrix, if not provided is the sum over each data matrix
#' @param theta optional value for a known overdispersion parameter
#' @param thresh threshold value for considering presence/absence of an edge respect to an estimated coefficient matrix
#' @param nCores number of cores over which parallelize the computation
#' @param nstab number of stability iterations to establish correct edges once selected a lambda
#' @param freq frequency threshold for the stability of the edges in real numbers between 1 and 100, by default is 85%
#' @param top quantile threshold for the top stable edges to retain in output matrix, between 1 and 100, default to 5
#' @param FPR threshold for the Meinshausen and Buhlmann formula on stability
#' @param verbose parameter controlling the output of the non parallelized part of the code
#' @param seed optional integer seed for reproducibility of the stability subsampling (only used when \code{lambda} is not a grid).
#'
#' @return A list whose structure depends on \code{lambda}:
#' \describe{
#'   \item{If \code{length(lambda) == 1}}{
#'     \itemize{
#'       \item \code{thresholded_network}: frequency-thresholded adjacency matrix
#'       \item \code{quantile_network}: top-quantile adjacency matrix
#'       \item \code{freq_matrix}: edge selection frequencies
#'       \item \code{lambda}: selected penalization value
#'       \item \code{call}: matched function call
#'       \item \code{freq_threshold_used}: frequency cutoff used
#'       \item \code{top_percent_used}: top-percent cutoff used
#'       \item \code{piMB}: Meinshausen–Bühlmann bound
#'     }
#'   }
#'   \item{If \code{length(lambda) > 1}}{
#'     \itemize{
#'       \item \code{network}: list of p x p sparse adjacency matrices (one per lambda)
#'       \item \code{coef_network}: array of unprocessed coefficients
#'       \item \code{lambda}: vector of penalization values
#'       \item \code{call}: matched function call
#'       \item \code{itermat}: EM iteration counts
#'       \item \code{warnings}: fitting warnings
#'     }
#'   }
#' }
#' @importFrom stats quantile
#' @importFrom Matrix forceSymmetric
#' @export
zinb_LGM_grp <- function(
    Xlist,
    lambda = NULL,
    conperc = 10,
    sym = c("AND", "OR"),
    offset = NULL,
    theta = NULL,
    thresh = 1e-6,
    nCores = 1,
    nstab = 50,
    freq = 85,
    top = 5,
    FPR = 0.05,
    seed = NULL,
    verbose = 0
) {
  # Match the symmetry argument
  sym <- match.arg(sym)
  fun_call <- match.call()

  # Checks on inputs:

  if (!is.list(Xlist) || length(Xlist) < 1) {
    stop("Xlist must be a non-empty list of matrices.")
  }

  if (!all(sapply(Xlist, function(x) is.numeric(x) && is.matrix(x)))) {
    stop("Each element in Xlist must be a numeric matrix.")
  }

  if (any(sapply(Xlist, ncol) < 2)){
    stop("Each matrix in Xlist must have at least 2 columns (variables).")
  }

  nvar <- ncol(Xlist[[1]])

  if (!all(sapply(Xlist, function(m) ncol(m) == nvar))) {
    stop("All matrices in Xlist must have the same number of columns.")
  }

  if (!is.null(offset) && (!is.numeric(offset) || length(offset) != sum(sapply(Xlist, nrow)))) {
    stop("Provided offset must be a numeric vector of the same length as the total number of rows in all Xlist matrices.")
  }

  if (!is.numeric(conperc) || conperc <= 1 || conperc > 100) {
    stop("conperc must be a number between 1 and 100.")
  }

  if (!is.numeric(freq) || freq <= 0 || freq > 100) {
    stop("freq must be a number between 0 and 100.")
  }

  if (!is.numeric(top) || top <= 0 || top > 100) {
    stop("top must be a number between 0 and 100.")
  }

  if (!is.null(lambda) && any(lambda < 0)) {
    stop("lambda must be non-negative values")
  }

  # Seed for reproducibility (used in stability sampling below)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) stop("seed must be a single numeric/integer value.")
    set.seed(as.integer(seed))
  }

  # Initialize offset if null
  if (is.null(offset)) {
    offset <- unlist(lapply(Xlist, function(i) log(rowSums(i))))
    if (verbose > 1) cat("Offset computed from row sums\n")
  }

  # Create indexes for matrices
  mat_indexes <- indexes(Xlist, nvar)
  ratio <- ifelse(any(sapply(Xlist, function(i) dim(i)[1] < nvar)), 0.05, 0.0001)
  penalty <- "LASSO"


  # Lambda computation

  if (is.null(lambda)) {
    if (verbose > 0) cat("\nSearching lambda...\n")
    # Maximum possible number of undirected edges in a fully connected graph
    maxcon <- nvar * (nvar - 1) / 2
    lambda <- find_lammax_grp(Xlist)

    # Initialize cycle control parameter
    ctl = 0
    # Find lambda to stabilize
    repeat {
      res <- adjust_lambda(
        current_lambda = lambda,
        Xlist = Xlist,
        mat_indexes = mat_indexes,
        nvar = nvar,
        sym = sym,
        theta = theta,
        thresh = thresh,
        nCores = nCores,
        offset = offset
      )

      # Setting threshold for the bisection process
      lower_thresh <- max(1, (conperc - 2) * (maxcon / 100))
      upper_thresh <- min(maxcon, (conperc + 2) * (maxcon / 100))

      if (res$edges >= lower_thresh && res$edges <= upper_thresh) break  # stop in desired range

      if (res$edges < lower_thresh) {
        lambda <- lambda / 2 # halve lambda_max to have at least (conperc-2)% of total edges
      }

      if (res$edges > upper_thresh) {
        lambda <- lambda * 2 # double lambda to have less than (conperc +2)% of total edges
        break # exit the loop as we're probably close to the desired connectivity for lambda
      }

      if (ctl > 100) {
        warning("Lambda Cycle hitted maximum iterations. Using now the last lambda")
        break
      }

      if (lambda < 1e-10 || lambda > 1e10) {
        stop("Lambda diverged outside of safe range during connectivity tuning. Check data sparsity or conperc setting.")
      }

      ctl <- ctl + 1

    }

      if (verbose > 0) cat("\nLambda search complete\n")

  }

  # We should have a lambda now if not provided, therefore:
  nlambda = length(lambda)


  # When given a lambda grid calculates results over the grid
  if (nlambda > 1) {

    cat(
      "Learning for NBII graphical model\n",
      "nlambda: ", nlambda, "\n",
      "penalty function: ", penalty, "\n",
      sep = ""
    )
    net <- zinb_LGM_net_grp(
      X = mat_indexes$matricione,
      lambda = lambda,
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
      verbose = verbose
    )

    if (verbose > 0) cat("Full process complete\n\n")

    # Output results
    list(
      network = net$hat_net,
      coef_network = net$coef_net,
      lambda = lambda,
      call = fun_call,
      itermat = net$itermat,
      warnings = net$warnings
    )
  } else {
    # When a grid of lambda is not given, we exploit the fact that we calculated a best lambda
    # in a desired range of connectivity, now we perform stability selection on that lambda
    cat("\nStabilizing ",
        nstab,
        " times\n", "for lambda = ",
        lambda,
        "\n",
        sep = ""
        )
    # let's initialise an object to sum over the nstab output matrices
    summat <- Matrix::Matrix(0, nvar, nvar, sparse = TRUE)
    #let's initialise also a vector that will accomodate the number of edges estimated
    # at each iteration and that we need for the
    #  meinshausen and buhlmann formula of false positive rate
    ql <- 0
    for(s in 1:nstab){
      # Extracting random indices corresponding to half of the samples of each data matrix
      sub_ind <- lapply(Xlist, function(i){sample(nrow(i), size = floor(nrow(i)/2))})
      # Create the two complementary sublists
      # (map is a function that allows to perform coordinated operations on two lists)
      sub_Xlist_1 <- Map(function(m, i) m[i, , drop = FALSE], Xlist, sub_ind)
      sub_Xlist_2 <- Map(function(m, i) m[-i, , drop = FALSE], Xlist, sub_ind)
      # let's use an helper function to run the method on the sublists
      run_half <- function(sub_Xlist) {
        offset <- unlist(lapply(sub_Xlist, function(i) log(rowSums(i) + 1)))
        sub_mat_indexes <- indexes(sub_Xlist, nvar)

        zinb_LGM_net_grp(
          X = sub_mat_indexes$matricione,
          lambda = lambda,
          sym = sym,
          theta = theta,
          thresh = thresh,
          nCores = nCores,
          p = nvar,
          totrow = sub_mat_indexes$totrow,
          totcol = sub_mat_indexes$totcol,
          P = sub_mat_indexes$P,
          k = length(sub_Xlist),
          nk = sub_mat_indexes$nk,
          inds = sub_mat_indexes$inds,
          group = sub_mat_indexes$group,
          ind = sub_mat_indexes$ind,
          ind_j = sub_mat_indexes$ind_j,
          offset = offset,
          verbose = verbose
        )
      }

      res1 <- run_half(sub_Xlist_1)$hat_net[[1]]
      res2 <- run_half(sub_Xlist_2)$hat_net[[1]]

      summat <- summat + res1 + res2

      res1m <- as.matrix(res1)
      res2m <- as.matrix(res2)

      ql <- ql +
        sum(upper.tri(res1m) & (res1m != 0)) +
        sum(upper.tri(res2m) & (res2m != 0))

      if (verbose > 0) cat("\nStability run number ", s, " performed\n", sep = "")
    }
    # Let's compute the final ql
    ql <- ql / (2 * nstab)

    # Now we need to extract the final network
    # Firstly according to frequence
    freq_mat <- summat/(2*nstab)
    freq <- if (!exists("freq")) 0.85 else freq/100

    # Create thresholded network in-place from freq_mat
    threshold_net <- freq_mat
    threshold_net@x <- as.numeric(threshold_net@x >= freq)


    # Now according to top x% most stable edges
    top <- if (!exists("top")) 0.05 else top / 100

    # Extract upper triangle values (without sorting full matrix)
    upper_tri_indices <- which(upper.tri(freq_mat), arr.ind = TRUE)
    upper_vals <- freq_mat[upper_tri_indices]

    # And now we search the cutoff for top x% most stable edges
    n_top_edges <- ceiling(length(upper_vals) * top)
    if (n_top_edges < 1) n_top_edges <- 1  # handle edge case
    quantile_cutoff <- quantile(upper_vals, probs = 1 - top)

    # Here we build quantile-based binary matrix efficiently
    quantile_net <- freq_mat
    quantile_net@x <- as.numeric(quantile_net@x >= quantile_cutoff)

    # Symmetrization in-place (for large sparse matrices, avoid full + transpose)
    threshold_net <- Matrix::forceSymmetric(threshold_net, uplo = "U")
    quantile_net <- Matrix::forceSymmetric(quantile_net, uplo = "U")
    diag(threshold_net) <- 0
    diag(quantile_net) <- 0

    # Meinshausen and Buhlmann formula:

    piMB <- 0.5*(1 + (1/FPR)*((ql/(nvar*(nvar-1)/2))^2))

    if (verbose > 0) cat("Stability selection complete\n\n")
    list(
      thresholded_network = threshold_net,
      quantile_network = quantile_net,
      freq_matrix = freq_mat,
      lambda = lambda,
      call = fun_call,
      freq_threshold_used = freq,
      top_percent_used = top,
      piMB = piMB
    )
  }

}
#' Fit local grouped graphical-model regressions parallelized across variables (internal)
#'
#' Runs the local (nodewise) grouped model estimation for each variable \code{j = 1,...,p},
#' parallelizing over \code{j}. For each \code{j}, a wrapper is called that fits the model
#' across the provided \code{lambda} grid and returns coefficients and diagnostics.
#' The resulting grouped coefficients are then aggregated into \code{p x p} joint
#' coefficient matrices (one per lambda) and thresholded into adjacency matrices.
#'
#' @param X Design matrix (typically sparse), constructed by \code{indexes()}.
#' @param lambda Numeric vector of penalization values.
#' @param sym Character; symmetrization rule for adjacency matrices, one of
#'   \code{"AND"} or \code{"OR"}.
#' @param theta Optional numeric overdispersion parameter passed downstream.
#' @param thresh Numeric threshold used to binarize coefficient matrices into edges.
#' @param nCores Integer; number of cores to use for parallelization over variables \code{j}.
#' @param p Integer; number of variables/genes.
#' @param totrow Integer; total number of rows in \code{X}.
#' @param totcol Integer; total number of columns in \code{X} excluding intercept columns
#'   (as constructed in \code{indexes()}).
#' @param P Integer; total number of columns in the full design matrix including intercepts.
#' @param k Integer; number of datasets / matrices integrated.
#' @param nk Integer vector of length \code{k}; number of samples (rows) in each dataset.
#' @param inds List of length \code{k}; row indices mapping each dataset into the stacked design.
#' @param group Integer vector defining group memberships for grouped-lasso style penalties.
#' @param ind Integer vector indexing columns to reorder the design according to group structure.
#' @param ind_j Integer vector used to map grouped coefficients back into \code{k x p} blocks.
#' @param offset Optional numeric vector of offsets (length \code{totrow}).
#' @param verbose Integer verbosity level.
#'
#' @return A list with:
#' \describe{
#'   \item{coef_net}{Array of raw grouped coefficients with dimensions \code{(k*p) x p x nlambda}.}
#'   \item{hat_net}{List of length \code{nlambda} of sparse \code{p x p} adjacency matrices.}
#'   \item{warnings}{List of length \code{p} storing warning metadata per fitted variable.}
#'   \item{itermat}{List of length \code{p} containing iteration counts across \code{lambda}.}
#' }
#'
#' @keywords internal
#' @noRd
zinb_LGM_net_grp = function(X,
                            lambda = NULL,
                            sym = c("AND", "OR"),
                            theta = NULL,
                            thresh = 1e-6,
                            nCores = 1,
                            p,
                            totrow,
                            totcol,
                            P,
                            k,
                            nk,
                            inds,
                            group,
                            ind,
                            ind_j,
                            offset = NULL,
                            verbose = 0)
  {

  nlambda = length(lambda) # salviamo la dimensione della griglia

  # la seguente matrice coef_mat deve accomodare i coefficienti risultanti escluse le intercette quindi saranno da ciascun fit k*p, per ogni gene quindi per p, per il numero di lambda nlambda
  coef_mat = array(dim = c(k*p, p, nlambda))

  # creiamo anche una lista ove possiamo salvare gli eventuali warning che si possono aver avuto dai fit con il corrispettivo valore di j

  warnings <- vector("list", p)

  # Inizializziamo una matrice di p colonne che accomodi i risultati del numero di iterazioni ottenute per ogni lambda, questa ci serve come controllo per capire quante iterazioni fa il nostro metodo

  itermat <- vector("list", p) # matrix(nrow = nlambda, ncol = p)

    # Ora lanciamo la parallelizzazione
  if (!is.numeric(nCores) || nCores < 1) {
    stop("nCores must be a positive integer.")
  }


  # inizializziamo i cores
  cl <- parallel::makeCluster(nCores)
  on.exit(parallel::stopCluster(cl))

  # Variabili da esportare
  vars_to_export <- c(
    "X", "lambda", "theta", "thresh",
    "p", "k", "P", "nk", "inds", "ind", "ind_j", "group", "nlambda",
    "verbose",
    "zinb_LGM_wrapper_grp", "dNBI", "dNBII",
    "grp_nb2_objective", "sigma_ml",
    "zilgm_negbin2_grp", "hat_net", "per_lambda_normal", "grpnet.default"
  )

  # esportiamo variabili necessarie al ciclo
  parallel::clusterExport(cl, varlist = vars_to_export, envir = environment())

  # carichiamo librerie necessarie
  parallel::clusterEvalQ(cl, {
    library(Matrix)
    library(igraph)
    library(MASS)
    library(parallel)
  })
    # codice per debug
  # coef_tmp <- lapply(1:p, function(j) {
  #   jth <- j + seq(k, k * p, by = p)
  #   ind_1 <- ind[-which(ind == jth)]
  #   group_j <- group[-which(group == j + k)]
  #
  #   zinb_LGM_wrapper_grp(
  #     j = j,
  #     X = X,
  #     jth = jth,
  #     lambda = lambda,
  #     theta = theta,
  #     nk = nk,
  #     inds = inds,
  #     group = group_j,
  #     k = k,
  #     P = P,
  #     ind = ind_1,
  #     ind_j = ind_j,
  #     offset = offset,
  #     p = p,
  #     nlambda = nlambda,
  #     thresh = thresh,
  #     verbose = verbose
  #   )
  # })

  # passiamo esplicitamente gli argomenti a cluster apply
  coef_tmp = parallel::clusterApply(cl = cl, 1:p, function(j,
                                                           group,
                                                           X,
                                                           lambda,
                                                           offset,
                                                           theta,
                                                           nk,
                                                           inds,
                                                           k,
                                                           P,
                                                           ind,
                                                           ind_j,
                                                           p,
                                                           nlambda,
                                                           thresh,
                                                           verbose) {
    # Prendiamo gli indici di Y
    jth = j + seq(k, k*p, by = p) # indice j di ogni matrice nel matricione
    ind_1 = ind[-which(ind == jth)] # indice di tutti tranne i jth su cui fare regressione
    group_j = group[-which(group == j+k)] # stesso discorso togliamo i jth dai gruppi
    # lanciamo la seconda shell del nostro metodo: zilgm_wrapper_grp che a j fissato cicla su lambda. Questa funzione contiene anche l'implementazione di active set

    # NB penalty.factor e NULL, lo useremo per introdurre il termine della penalty di gruppo sqrt(k) piu avanti ma va considerato variabile interna
    zinb_LGM_wrapper_grp(j = j,
                         X = X,
                         jth = jth,
                         lambda = lambda,
                         theta = theta,
                         nk = nk,
                         inds = inds,
                         group = group_j,
                         k = k,
                         P = P,
                         ind = ind_1,
                         ind_j = ind_j,
                         offset = offset,
                         p = p,
                         nlambda = nlambda,
                         thresh = thresh,
                         verbose = verbose)
  },          group = group,
  X = X,
  lambda = lambda,
  theta = theta,
  nk = nk,
  inds = inds,
  k = k,
  P = P,
  ind = ind,
  ind_j = ind_j,
  p = p,
  nlambda = nlambda,
  offset = offset,
  thresh = thresh,
  verbose = verbose)


  # Qui andiamo a popolare la matrice dei coefficienti e quelle dei warnings e iterazioni
  for (j in 1:p) {
    itermat[[j]] <- coef_tmp[[j]]$niterations
    coef_mat[, j, ] = as.matrix(coef_tmp[[j]]$Bmat)
    warnings[[j]] <- j
    if(!is.null(coef_tmp[[j]]$warning_lambda)){

      warnings[[j]] <- c(warnings[[j]],coef_tmp[[j]]$warning_lambda)
    }
  }

  # NB coef_mat non e come in Park, direttamente la matrice (p)x(p), e (k*p)x(p) e va trasformato
  # usiamo questo indice di gruppo
  indgrp <- split(ind_j, group[-(1:k)])

  # Conserveremo i risultati nella seguente matrice:
  coef_fin <- array(dim = c(p, p, nlambda))

  # riempiamo la matrice coef_fin
  for (l in 1:nlambda) {
    # Compute and store the matrix for the current lambda
    coef_fin[,,l] <- per_lambda_normal(coef_mat, indgrp, p, l)
  }


  # Futuro codice per costruire le reti:
  ghat = lapply(1:nlambda, FUN = function(l) hat_net(coef_fin[, , l], thresh = thresh, type = sym))
  gs = lapply(1:nlambda, FUN = function(l) Matrix::Matrix(ghat[[l]]))

  # componiamo la lista risultato
  result_list <- list(coef_net = coef_mat, hat_net = gs, warnings = warnings, itermat = itermat)
  #return(list(hat_net = gs, coef_net = coef_mat)) # per quando avremo direttamente la rete in uscita, per ora i coefficienti
  return(result_list)
}
#' Fit one nodewise model across a lambda grid (internal)
#'
#' For a fixed target variable \code{j}, this helper fits the local grouped model
#' across the provided \code{lambda} grid and stores (i) intercepts and
#' (ii) grouped coefficients (excluding intercepts) in sparse matrices.
#'
#' Warning and error conditions are recorded by lambda; on error, the corresponding
#' column is left as zeros and the lambda is recorded.
#'
#' @param j Integer; target node index in \code{1:p}.
#' @param X Design matrix.
#' @param jth Integer vector of length \code{k}; column indices in \code{X} corresponding
#'   to the target response for each dataset/block.
#' @param lambda Numeric vector of penalization values.
#' @param theta Optional numeric overdispersion parameter passed downstream.
#' @param nk Integer vector of length \code{k}; sample sizes per dataset.
#' @param inds List of length \code{k}; row indices per dataset in the stacked design.
#' @param group Integer vector of group memberships for the predictors (excluding the target).
#' @param k Integer; number of datasets / blocks.
#' @param P Integer; total number of columns in the full design matrix including intercepts.
#' @param ind Integer vector of predictor column indices used in the regression for this node.
#' @param ind_j Integer vector mapping grouped coefficients back into the full \code{(P-k)} layout.
#' @param offset Optional numeric vector of offsets.
#' @param p Integer; number of variables/genes.
#' @param nlambda Integer; length of \code{lambda}.
#' @param thresh Numeric threshold forwarded to the solver.
#' @param verbose Integer verbosity level; if \code{>0}, prints progress.
#'
#' @return A list with:
#' \describe{
#'   \item{b0}{Sparse \code{k x nlambda} matrix of intercepts.}
#'   \item{Bmat}{Sparse \code{(P-k) x nlambda} matrix of coefficients (excluding intercepts).}
#'   \item{warning_lambda}{Numeric vector of lambda values that triggered warnings and/or errors.}
#'   \item{niterations}{Numeric vector of length \code{nlambda} with EM iteration counts (may be shorter if early-stopped).}
#' }
#'
#' @keywords internal
#' @noRd
zinb_LGM_wrapper_grp = function(j, X,
                                jth,
                                lambda, theta,
                                nk, inds, group,
                                k, P, ind, ind_j, #aggiunto
                                offset,
                                p, nlambda, thresh,
                                verbose)
  {
  # j e FISSATO, CICLIAMO SU LAMBDA

  # la seguente matrice Bmat deve accogliere i coefficienti del fit ad ogni iterazione j e per ogni lambda.

  # la matrice Bmat quindi serve per accogliere i coefficienti ottenuti dal fit, ma anche quelli nulli diagonali che non si ottengono dal fit, i quali devono rimanere infatti nulli
  # b0 invece accogliera le intercette

  # e.g. per p = 500, dal fit, quindi escludendo le intercette, ne otteniamo 1497 a gruppi i.e. 123 = Interc, 456 = 1° gruppo, 789 2° gruppo etc...

  # Bmat ne deve sempre accomodare 1500

  # Bmat = Matrix(0, p, nlambda, sparse = TRUE)
  Bmat = Matrix::Matrix(0, P - k, nlambda, sparse = TRUE) # Matrix matrice sparsa

  # b0 deve sempre accomodare i primi k (le intercette)
  b0 = Matrix::Matrix(0, k, nlambda, sparse = TRUE)

  # usiamo il seguente vettore per salvare il beta stimato dall'EM a ciascun lambda e passarlo come beta iniziale nel ciclo ai lambda successivi

  betaprev = NULL

  # Inizializziamo il vettore atto a raccogliere i warning
  warning_lambda <- c()

  # Inizializziamo anche il vettore preposto a raccogliere il numero di iterazioni per un dato lambda, fissato j
  niterations <- c()

  if (length(ind) == 0) { # Se non ci sono colonne selezionate restituisci risultato vuoto
    Bmat = Bmat
    b0 = b0
  } else {
    if(k == 1){
      y = X[,jth]
    } else {
      y = Matrix::rowSums(X[,jth])
    }
    for (iter in 1:nlambda) {
      # Primo lambda: calcoliamo tutto agnosticamente
      if (iter == 1){
        # printiamo sempre if (verbose == 1) {
        cat("lambda = ", lambda[iter], ", ", j, "/", p, "th node learning \n", sep = "")
        # }

        # tryCatch per gestire i lambda problematici, ossia dove eventualmente il fit da errori
        coef_res <- tryCatch({
          fitem = zilgm_negbin2_grp(y = y, x = X[, ind, drop = FALSE],
                                    k,
                                    nk, inds, group,
                                    lambda = lambda[iter], theta = theta,
                                    bvec0 = NULL,
                                    offset = offset,
                                    thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
        }, warning = function(w) {
          # salviamo quale lambda sta eventualmente dando problemi per questo j
          cat("Warning NAN nel fit su lambda:", lambda[iter], "\n")

          # Restituiamo comunque il risultato del fit
          fitem <- zilgm_negbin2_grp(y = y, x = X[, ind, drop = FALSE],
                                     k,
                                     nk, inds, group,
                                     lambda = lambda[iter], theta = theta,
                                     bvec0 = NULL,
                                     offset = offset,
                                     thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
          fitem$warning_lambda = lambda[iter]
          return(fitem)
        })
        if(!is.null(coef_res$warning_lambda)){
          warning_lambda = c(warning_lambda, coef_res$warning_lambda)
        }

        # riempiamo solo con i coefficienti calcolati, lasciamo nulli gli altri
        Bmat[ind_j[ind_j != (jth-k)],iter] = coef_res$bvec[-(1:k)]

        # Aggiungiamo una condizione di stop
        if (sum(Bmat[,iter] == 0) == length(Bmat[,iter])){
          break}

        # Gli active sono proprio quelli diversi da zero in Bmat,dobbiamo pero "tradurli" nell'ordine del fit di y VS X
        # Non usiamo Active, o almeno manteniamo la taglia del sistema
        # Active <- as.numeric(c(1:k, sapply(which(Bmat[1:p, iter] != 0), function(i){
          #k + seq(from = i, to = k * p, by = p)
        #})))

        # Non aggiorniamo dunque i gruppi del prossimo fit: in quanto vogliamo
        #  conservare la taglia del problema
        # group = c(1:k, rep((k+1):(k+length(which(Bmat[1:p, iter] != 0))), each = k))

        # Salviamo il beta nuovo da dare in input al prossimo step:
        # NO! betaprev = c(coef_res$bvec[1:k], coef_res$bvec[-(1:k)][coef_res$bvec[-(1:k)]!=0])
        # In questo caso vogliamo conservare lo stesso identico beta ottenuto dal fit

        # Ulteriore modifica : non salviamo alcun betaprev, vogliamo che sia sempre vuoto
        # la dimensionalita deve restare invariata
        # betaprev = coef_res$bvec

        # Riempiamo anche il vettore delle intercette
        b0[,iter] = coef_res$bvec[1:k]

        # E salviamo quante iterazioni sono state necessarie
        niterations[iter] <- coef_res$iterations

      } else { # le altre iterazioni partono con un beta noto
        # printiamo sempre if (verbose == 1) {
        cat("lambda = ", lambda[iter], ", ", j, "/", p, "th node learning \n", sep = "")
        # }

        # tryCatch per gestire i lambda problematici
        coef_res <- tryCatch({
          fitem = zilgm_negbin2_grp(y = y, x = X[, ind, drop = FALSE], # ind e cambiato con gli active
                                    k,
                                    nk, inds, group,
                                    lambda = lambda[iter], theta = theta,
                                    bvec0 = NULL,
                                    offset = offset,
                                    thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
        }, warning = function(w) {
          # salviamo quale lambda sta dando problemi per questo j
          cat("Warning NAN nel fit su lambda:", lambda[iter], "\n")

          # Restituiamo comunque fitem
          fitem <- zilgm_negbin2_grp(y = y, x = X[, ind, drop = FALSE], # ind e cambiato con gli active
                                     k,
                                     nk, inds, group,
                                     lambda = lambda[iter], theta = theta,
                                     bvec0 = NULL,
                                     offset = offset,
                                     thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
          fitem$warning_lambda = lambda[iter]
          return(fitem)
        })
        if(!is.null(coef_res$warning_lambda)){
          warning_lambda = c(warning_lambda, coef_res$warning_lambda)
        }

        # Qui dobbiamo replicare il paradigma di active set

        # Non useremo aggiornare solo gli active ma tutti i beta dunque:
        #non Bmat[(Active[-c(1:k)] - k), iter] = coef_res$bvec[-(1:k)], ma:
        Bmat[ind_j[ind_j != (jth-k)],iter] = coef_res$bvec[-(1:k)]

        # Aggiungiamo una condizione di stop
        if (sum(Bmat[,iter] == 0) == length(Bmat[,iter])){
          break}

        # Di nuovo non consideriamo Active
        # Active <- as.numeric(c(1:k, sapply(which(Bmat[1:p, iter] != 0), function(i){
        #   k + seq(from = i, to = k * p, by = p)
        # })))

        # Nuovamente non aggiorniamo dunque i gruppi del prossimo fit:
        # group = c(1:k, rep((k+1):(k+length(which(Bmat[1:p, iter] != 0))), each = k))

        # Inoltre betaprev deve rimanere il risultato del fit
        # betaprev = c(coef_res$bvec[1:k], coef_res$bvec[-(1:k)][coef_res$bvec[-(1:k)]!=0]) # check intercette non nulle

        # anche qui annulliamo il betaprev
        # betaprev = coef_res$bvec


        b0[,iter] = coef_res$bvec[1:k]
        niterations[iter] <- coef_res$iterations
      }
    }
  }
  return(list(b0 = b0, Bmat = Bmat, warning_lambda = warning_lambda, niterations = niterations))
}
#' Inner EM solver for the grouped NB2 ZI graphical model (internal)
#'
#' Fits model parameters for a fixed response vector \code{y} and predictor matrix \code{x}
#' at a single penalization value \code{lambda}. When \code{y} contains zeros, an EM
#' procedure alternates between:
#' \itemize{
#'   \item E-step: updating latent zero-inflation responsibilities \code{z} and per-dataset
#'     zero-inflation probabilities \code{prob}
#'   \item M-step: fitting a penalized Poisson grouped regression via \code{grpnet.default}
#'     (weights modified by \code{1 - z}) and updating the NB2 overdispersion parameter
#' }
#'
#' When \code{y} has no zeros, the function fits only the penalized Poisson regression and
#' estimates overdispersion directly (unless \code{theta} is fixed).
#'
#' @param y Numeric response vector of counts.
#' @param x Numeric design matrix (predictors), excluding the target columns.
#' @param k Integer; number of datasets/blocks.
#' @param nk Integer vector of length \code{k}; sample sizes per dataset.
#' @param inds List of length \code{k}; row indices for each dataset in \code{y/x}.
#' @param group Integer vector of group labels for grouped penalties.
#' @param theta Optional numeric scalar; if provided, the NB2 overdispersion is held fixed.
#' @param lambda Numeric scalar penalization value.
#' @param bvec0 Optional numeric initial coefficient vector for the Poisson fit.
#' @param offset Optional numeric vector of offsets (length \code{NROW(x)}).
#' @param thresh Numeric threshold parameter carried through the call chain (not used to zero coefficients in the current implementation).
#' @param EM_tol Numeric convergence tolerance for the EM objective.
#' @param EM_iter Integer maximum number of EM iterations.
#' @param tol Numeric tolerance passed to the grouped regression solver.
#' @param maxit Integer maximum iterations for the grouped regression solver.
#' @param verbose Integer verbosity level (used only for selected messages in the current implementation).
#'
#' @return A list of class \code{"zilgm"} with components:
#' \describe{
#'   \item{lambda}{Penalization value used.}
#'   \item{bvec}{Estimated coefficient vector from the penalized Poisson regression.}
#'   \item{theta}{Estimated (or fixed) NB2 overdispersion parameter.}
#'   \item{prob}{Estimated per-dataset zero-inflation probabilities.}
#'   \item{pos_zero}{Indices of zero observations in \code{y}.}
#'   \item{iterations}{Number of EM iterations performed (0 in the no-zero case).}
#'   \item{loglik}{Final objective value used for EM convergence monitoring.}
#'   \item{call}{The matched call.}
#'   \item{nzeroel}{Number of nonzero coefficients reported by the regression solver.}
#' }
#'
#' @keywords internal
#' @noRd
zilgm_negbin2_grp = function(y, x, k, # dati e dimensioni
                             nk, inds, group, # indici e gruppi
                             theta = NULL, # dispersione
                             lambda, # penalty sparsita
                             bvec0 = NULL, # betahat del precedente lambda
                             offset = NULL,
                             # parametri tolleranza ciclo EM:
                             thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
  {
  fun_call = match.call() # salva la chiamata alla funzioneper mostrarla in output

  # salviamo alcune variabili utili
  k = k
  N = NROW(x)
  P = NCOL(x) # 1500 quando privata delle jth-esime colonne
  p = (P - k )/k # era 499 quando chiamata con p = p-1

  if (!is.null(theta)) {
    fixed_theta = TRUE
    init_theta = theta
  } else {
    fixed_theta = FALSE
  }

  # Indici di elementi uguali a zero in y, necessario per fare inferenza su provenienza dalla delta o da NB
  pos_zero <- lapply(inds, function(i) y[i] == 0)
  pos_nzero <- lapply(pos_zero, function(x) !x)
  # Inizializzazione variabili latenti
  z <- lapply(nk, function(x) rep(1e-6, x))

  # Questa variabile la usiamo per la penalty lasso, sara la radice della cardinalita dei gruppi che moltiplica lambda, NB le intercette ne sono prive

  penalty.factor <- c(rep(0, k), rep(sqrt(k), p))


  w = rep(1, N)

  for(i in 1:k){ # deve essere 1/n_h per ogni h-esimo esperimento
    w[inds[[i]]]= w[inds[[i]]]/length(inds[[i]])
  }# fattore di scala, si puo elidere ed lasciare che sia "incorporato" su lambda se costante

  # Qui di seguito verifichiamo se i valori di y sono tutti uguali
  # In tal caso, cioe unique(y) == 1, il valore costante non potrebbe essere zero altrimenti avremmo filtrato quel gene in uno step di preprocessing
  # quindi qui si tiene solo conto del caso costante > 0 che comporterebbe prob = 0
  # come ulteriore conseguenza necessaria pos_zero sara tutto falso

  if (length(unique(y)) == 1) {
    pos_zero <- lapply(pos_zero, function(x) rep(FALSE, length(x)))
    param = list(bvec = rep(0, P), sigma = 0, # bvec deve uscire 1500
                 # andrebbero inserite anche le intercette diverse da 0
                 prob = rep(0, k), pos_zero = pos_zero,
                 iterazio = 0)
    return(param)
  }

  # Adesso passiamo alla inizializzazione di parametri relativi al fit glm
  # Una stima iniziale del vettore delle medie:
  mu0 <- lapply(inds, function(i) {
    mu_ind <- y[i][y[i] > 0]
    rep(mean(mu_ind), length(i))
  })

  # Useremo beta in modo iterativo, innanzitutto potrebbe venir passato dalla stima del lambda precedente, qualora non fosse cosi pero va inizializzato:
  if (is.null(bvec0)){
    bvec0 <- rep(0,(p+1)*k)
  } else {
    bvec0 = bvec0
  }
  # Passiamo adesso ad inizializzare il parametro che tiene conto della dispersione
  # Qui lo chiameremo theta come nella nostra formalizzazionema si riferisce in realta all'alpha dell'articolo di park cioe il rapporto tra la media mu e il parametro di dispersione theta
  # con un theta alto questo rapporto deve risultare piccolo, la condizone iniziale e dunque simil poissoniana
  theta0 = 1e-4

  # Inizializziamo infine delle probabilita di avere zeri provenienti dalle delta
  prob0 <- sapply(1:k, function(i) {
    temp_prob <- (sum(pos_zero[[i]]) - sum(dNBII(0, mu = mu0[[i]], sigma = theta0, log = FALSE))) / nk[i]
    pmin(1, pmax(1e-10, temp_prob))
  })

  # Massimo della funzione obiettivo per la convergenza del ciclo EM
  erisk_prev = 1e+150


  iterazio = 0
  # Questo if apre al caso in cui nel vettore y non ci siano zeri, in tal caso non facciamo alcuna inferenza su prob, stimiamo direttamente media e dispersione
  if (sum(sapply(pos_zero, sum) == 0)) {
    ######################## CASO DA ADATTARE
    print("no zeroes in Y")
    # fit poisson per stimare beta
    sol_bvec = grpnet.default(x = x, y = y, group = group, weights = w,
                                           family = "poisson", beta = bvec0,
                                           intercept = FALSE, lambda = lambda,
                                           standardize = FALSE, penalty = "LASSO",
                                           penalty.factor = penalty.factor,
                                           offset = offset,
                                           thresh = tol, maxit = maxit)

    nzeroel = sol_bvec$nzcoef # elementi non nulli
    bvec = sol_bvec$beta
    mu = sol_bvec$mu
    # se il fit ha dato problemi di NaN nei coefficienti teniamo il beta precedente
    # Inoltre stimiamo theta con la formula inversa da media e varianza
    if (any(is.nan(bvec))){
      warning("NaN in the Poisson glm coefficients")
      bvec = bvec0
      theta <- (stats::var(y) - mean(y)) / (mean(y)^2)
      prob = prob0
      erisk = erisk_prev
      nzeroel = sum(bvec0 != 0)
      mu = exp(offset + (x%*%bvec)/w)
    } else {
      # Altrimenti calcoliamo la nuova media dai beta
      # e poi la nuova dispersione
      #mu = exp(offset + (x%*%bvec)/w)
      if (fixed_theta) {
        theta = init_theta
      } else {
        theta = sigma_ml(y = y, mu = mu, weights = w)
      }
      # restituiamo l'output
      prob = prob0
      iterazio = 0
      erisk = grp_nb2_objective(y = y, prob = prob0, bvec = bvec, mu = mu, lambda = lambda,
                                weights = w, penalty.factor = penalty.factor, sigma = theta,
                                posz = unlist(pos_zero), k = k, p = p, inds = inds)
      # erisk = erisk_prev
      # theta = theta
      # bvec0 = bvec
    }

    ########################################
    # DOPO L'ELSE IL CICLO EM DI INTERESSE #
    ########################################



  } else {
    # Qui entriamo nel vero e proprio EM dove stimiamo iterativamente coefficienti, dispersione e probabilita di delta
    for (iterazio in 1:EM_iter) {
      print(iterazio)
      # Inizializzazione E-step:
      prob = NULL
      tmp_z <- vector("list", k)

      # E-step: calcoliamo quanto sia probabile che uno zero venga dalla delta data l'attuale media mu0, il parametro theta0 (sempre rapporto tra mu0 e theta0 della NBI I vedi commento inizializzazione di theta0) e data l'attuale stima della prob
      # Il seguente codice sarebbe l'equazione 6 quando si usa pi invece che exp(fi)
      for (i in 1:k) {
        tmp_z[[i]] <- prob0[i] / (prob0[i] + (1 - prob0[i]) * dNBII(0, sigma = theta0, mu = mu0[[i]], log = FALSE))
        tmp_z[[i]][is.nan(tmp_z[[i]])] <- 1
        tmp_z[[i]] <- pmin(1 - 1e-6, tmp_z[[i]])

        z[[i]][pos_zero[[i]]] <- tmp_z[[i]][pos_zero[[i]]]

        prob[i] <- sum(z[[i]]) / nk[i]
        prob[i] <- pmin(1, pmax(1e-10, prob[i]))
      }
      ###############################################
      print(paste("probability values:"))
      print(prob)
      ###############################################

      # M-step

      # Qui andiamo nello step di massimizzazione, partendo dal glm fit:
      # operiamo una grouped lasso poisson penalized regression
      # lo step precedente e tenuto in conto nei pesi w, che vengono aggiornati in modo tale da essere meno importanti nel fit quanto piu e probabile che vengano dalla delta in zero

      sol_bvec = grpnet.default(x = x, y = y, group = group, weights = w * (1 - unlist(z)),
                                             family = "poisson", beta = bvec0,
                                             intercept = FALSE, lambda = lambda,
                                             standardize = FALSE, penalty = "LASSO",
                                             penalty.factor = penalty.factor,
                                             offset = offset,
                                             thresh = tol, maxit = maxit)

      # Estraiamo i risultati del fit

      nzeroel = sol_bvec$nzcoef # elementi non nulli
      bvec = sol_bvec$beta # coefficienti stimati (il supporto sono i nodi)
      mu = sol_bvec$mu
      # In particolare qui abbiamo 499*3 = 1497 coefficienti e 3 intercette
      # mancano i coefficienti di j vs j
      if (any(is.nan(bvec))){
        warning("NaN in the Poisson glm coefficients")
        bvec = bvec0
        theta = (mean(y)^2)/(var(y)-mean(y))
        prob = prob0
        erisk = erisk_prev
        nzeroel = sum(bvec0 != 0)
        mu = exp(offset + (x%*%bvec)/(w * (1 - unlist(z))))
        break
      }

      #mu = exp(offset + (x%*%bvec)/(w * (1 - unlist(z)))) # nuovo vettore delle medie

      # Qui se il coefficiente alpha era noto lo si tiene fissato altrimenti viene stimato partendo dalle medie risultanti dal fit e dall'attuale z

      if (fixed_theta) {
        theta = init_theta
      } else {
        theta = sigma_ml(y = y, mu = mu, weights = w * (1 - unlist(z)))
      }

      # Calcoliamo la funzione obiettivo
      erisk = grp_nb2_objective(y = y, prob = prob, bvec = bvec, mu = mu, lambda = lambda,
                                weights = w, penalty.factor = penalty.factor, sigma = theta,
                                posz = unlist(pos_zero), k = k, p = p, inds = inds)

      # Controlli di convergenza:
      # Se e infinito usa il valore precedente
      if (is.infinite(erisk) | is.nan(erisk)) {erisk = erisk_prev}
      # nel caso in cui stiamo usando il valore precedente o non e cambiato entro una certa soglia fermiamo le iterazioni
      if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
        bvec = bvec
        theta = theta
        prob = prob
        z = z
        break
        # } else if (erisk > erisk_prev + 1e-10) {
        #   bvec = bvec0
        #   theta = theta
        #   prob = prob0
        #   break
      } else { # se la funzione obiettivo puo ancora crescere invece continua aggiornando i parametri:

        # erisk_prev = erisk
        # bvec0 = bvec
        # eta0 = eta
        # mu0 = mu
        # theta0 = theta
        # prob0 = prob

        erisk_prev = erisk
        mu0 = lapply(1:k, function(i){mu[inds[[i]]]})
        theta0 = theta
        prob0 = prob
        bvec0 = bvec
      }
    }
  }

  # flag = abs(bvec) < thresh
  # # stiamo annullando tutti i valori al di sotto di thresh!
  # bvec[flag] = 0
  #
  # raccogliamo i risultati da inserire nell'output
  out = list()
  out$lambda = lambda
  out$bvec = bvec
  out$theta = theta
  out$prob = prob
  out$pos_zero = which(unlist(pos_zero))
  out$iterations = iterazio
  out$loglik = erisk
  out$call = fun_call
  out$nzeroel = nzeroel
  class(out) = "zilgm"
  return(out)
}
