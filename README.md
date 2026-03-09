# JZILGM

**JZILGM** implements Joint Zero-Inflated Local Graphical Models on Negative Binomial distributed count data, with the purpose of data integration and network inference.

The method estimates association networks from count data variables using a **zero-inflated negative binomial graphical model** with **group penalties**, allowing joint inference across **multiple datasets**.

## Installation

Install directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("FranMigl/JZILGM")
```
Or install locally by downloading and then running in R:

```r
devtools::install("path/to/JZILGM")
```
---

## Input format

The main function expects a **list of count matrices**.

- rows = samples (can be different in each matrix)
- columns = variables (e.g. genes)
- all matrices must have the same number of columns

Example structure:

```r
ListMat <- list(
  matrix1,
  matrix2,
  matrix3
)
```
Each matrix typically represents **a different experiment or dataset generated from the same underlying network**.

---

## Example

```r
The following example illustrates joint network estimation from multiple matrices generated from the same ground-truth network.

library(JZILGM)
install.packages("netUtils")
library(netUtils)

# Generate the ground-truth network
net100 <- sample_lfr(n = 100, mu = 0.01, average_degree = 5,
                     max_degree = 10, min_community = 5, max_community = 20)

# Generate k = 3 data matrices from the same network
set.seed(1926)

p <- length(net100)
k <- 3

sim_data <- k_net_datagen(
  net   = net100,
  k     = k,
  seed  = 1926,
  n     = 250,
  p     = p,
  theta = 0.125,
  means = c(3, 6, 9),
  zerop = 0.4
)

# Check dimensions
lapply(sim_data, dim)

# Estimate the joint network
res <- zinb_LGM_grp(
  Xlist   = sim_data,
  sym     = "AND",
  nCores  = 1,
  verbose = 1
)

# Extract the inferred network
est_net <- res$thresholded_network

# Inspect the result
print(est_net)
```

---

## Output

The main function returns a list containing:

- `thresholded_network` — binary network using stability frequency threshold  
- `quantile_network` — binary network using top-percentage edge selection  
- `freq_matrix` — edge stability frequencies to be thresholded as you like 
- `lambda` — selected regularization parameter  
- `piMB` — stability selection bound from Meinshausen & Bühlmann  

---

## Citation

If you use this package, please cite:

Migliaccio F., Angelini C., Carissimo A., De Canditiis D. (2026).  
*scRNA-seq data integration via zero-inflated local negative binomial graphical model.*  
Manuscript in preparation.

### BibTeX

```
@article{migliaccio2026zilgm,
  title={scRNA-seq data integration via zero-inflated local negative binomial graphical model},
  author={Migliaccio, Francesco and Angelini, Claudia and Carissimo, Annamaria and De Canditiis, Daniela},
  year={2026},
  note={Manuscript in preparation}
}
```

---

## Status

This package is currently under active development.
