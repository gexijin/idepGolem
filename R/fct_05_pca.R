#' 05_pca
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd

## edit later
# Runs pathway analysis using PGSEA; this is copied and revised from PGSEA package
myPGSEA <- function(exprs, cl, range = c(25, 500), ref = NULL, center = TRUE,
                    p.value = 0.005, weighted = TRUE, nPermutation = 100, enforceRange = TRUE, ...) {
  if (is(exprs, "ExpressionSet")) {
    exprs <- exprs(exprs)
  }
  if (!is.list(cl)) {
    stop("cl need to be a list")
  }
  if (!is.null(ref)) {
    if (!is.numeric(ref)) {
      stop("column index's required")
    }
  }
  if (!is.null(ref)) {
    if (options()$verbose) {
      cat("Creating ratios...", "\n")
    }
    ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
    exprs <- sweep(exprs, 1, ref_mean, "-")
  }
  if (center) {
    exprs <- scale(exprs, scale = FALSE)
  } # column centering is done
  results <- matrix(NA, length(cl), ncol(exprs))
  rownames(results) <- names(cl)
  colnames(results) <- colnames(exprs)
  mode(results) <- "numeric"
  Setsize <- c(rep(0, length(cl))) # gene set size vector
  mean2 <- c(rep(0, length(cl))) # mean of the range of means
  meanSD <- c(rep(0, length(cl))) # SD of the range of means
  if (is.logical(p.value)) {
    p.results <- results
    mean.results <- results
  }
  for (i in 1:length(cl)) { # for each gene list
    # cat("\nProcessing gene set",i);
    if (class(cl[[i]]) == "smc") {
      clids <- cl[[i]]@ids
    } else if (class(cl[[i]]) %in% c("GeneColorSet", "GeneSet")) {
      clids <- cl[[i]]@geneIds
    } else {
      clids <- cl[[i]]
    }
    if (options()$verbose) {
      cat("Testing region ", i, "\n")
    }
    ix <- match(clids, rownames(exprs))
    ix <- unique(ix[!is.na(ix)])
    present <- sum(!is.na(ix))
    Setsize[i] <- present
    if (present < range[1]) {
      if (options()$verbose) {
        cat(
          "Skipping region ", i, " because too small-",
          present, ",\n"
        )
      }
      next
    }
    if (present > range[2]) {
      if (options()$verbose) {
        cat(
          "Skipping region ", i, " because too large-",
          present, "\n"
        )
      }
      next
    }
    texprs <- exprs[ix, ] # expression matrix for genes in gene set
    if (any(is.na(texprs))) {
      cat("Warning - 'NA' values within expression data, enrichment scores are estimates only.\n")
    }
    if (!is.matrix(texprs)) {
      texprs <- as.matrix(texprs)
    }

    stat <- try(apply(texprs, 2, t.test, ...))
    means <- try(apply(texprs, 2, mean, trim = 0.1)) # trim mean
    ps <- unlist(lapply(stat, function(x) x$p.value))
    stat <- unlist(lapply(stat, function(x) x$statistic))
    p.results[i, ] <- ps
    mean.results[i, ] <- means
    results[i, ] <- as.numeric(stat)

    # permutation of gene sets of the same size
    if (nPermutation > 2) { # no permutation if <=2
      meansRanges <- c(0, rep(nPermutation))
      for (k in 1:nPermutation) {
        ix <- sample.int(dim(exprs)[1], length(ix))
        texprs <- exprs[ix, ]
        means <- try(apply(texprs, 2, mean, trim = 0.1)) # trim mean
        meansRanges[k] <- dynamicRange(means)
      }
      mean2[i] <- mean(meansRanges)
      meanSD[i] <- sd(meansRanges, na.rm = TRUE) # NA are removed before calculating standard deviation
    }
  }
  return(list(results = results, p.results = p.results, means = mean.results, size = Setsize, mean2 = mean2, meanSD = meanSD))
}
