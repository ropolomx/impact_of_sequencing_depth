library(dplyr)
library(vegan)
amr_results_class <- amr_results %>% group_by(Sample, Class) %>% summarise(Hits=sum(Hits_Seen))


amr_results <- read.csv('AMR/amr_new_dataframe_ROP.csv')


amr_results_class_wide <- amr_results_class %>% spread(Class, Hits, fill = 0)
amr_results_class_wide <- amr_results_class %>% spread(Sample, Hits, fill = 0)
table(annotations$class)


amr_results_class_wide <- amr_results_class %>% spread(Class, Hits, fill = 0)

amr_results_class_mat <- amr_results_class_wide[,2:ncol(amr_results_class_wide)]
amr_class_norm <- data.frame(MRcounts(amr_class_experiment, norm=T))
amr_class_norm <- t(amr_class_norm)
row.names(amr_class_norm) <- str_replace(row.names(amr_class_norm), "X", "")


rarecurveROP <- function (x, step, sample, xlab = "Sample Size", ylab = "Species", 
          label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")?invisible
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)

  return(list(Nmax, Smax))
    
}


rarecurve_ROP <- function (x, step = 5, sample=raremax, label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  
  #plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
  #    type = "n", ...)
  #if (!missing(sample)) {
  #  abline(v = sample)
  #  rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
  #                                         y = z, xout = sample, rule = 1)$y)
  #  abline(h = rare, lwd = 0.5)
  #}
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    #lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  
  if (label) {
    samples <- as.character(row.names(x))
  }
  
  return(list(N=N, Out=out[[ln]], Samples=samples))
  #invisible(out)
}

rarefy_ggplot <- rarecurve_ROP(amr_results_class_mat)