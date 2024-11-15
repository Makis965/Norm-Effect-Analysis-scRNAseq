library(mclustcomp)
library(fpc)
library(cluster)
library(config)


# ---- statistical inference ----


conover_pair_comp <- function(current_gr, n, MS_within, r_mean) {
  # Function to calculate Conover test statistics for each pair comparison.
  current_n <- n[current_gr]
  SE <- sqrt(MS_within * sum(1/current_n))
  t <- diff(rev(r_mean[current_gr]))/SE
  return(c(t = t, SE = SE))
}

cohen_rank_d <- function(current_gr, n, r_mean, r_sd){
  
  n1 <- n[current_gr[1]]
  n2 <- n[current_gr[2]]
  n12 <- sum(n[current_gr])
  
  s1 <- r_sd[current_gr[1]]
  s2 <- r_sd[current_gr[2]]
  
  #pooled standard deviation standardized by the sum of observaton between samples
  sp <- sqrt(((n1-1)*(s1^2) + (n2-1)*(s2^2)) / (n12-2))
  
  #Cohen's d effect size for ranks 
  d <- (abs(r_mean[current_gr[1]] - r_mean[current_gr[2]])) / sp
  
  return(d)
}

conover_test <- function(values, groups, H, method = "bonferroni") {
  # Function for Conover's test of multiple comparisons using rank sums as post 
  # hoc test following a Kruskal-Wallis test.
  # H - Kruskal-Wallis test statistics
  
  # GROUP COMBINATIONS ---------------------------------------------------------
  gr <- as.character(groups)
  comb <- combn(unique(gr), 2) # all possible pairs of groups
  comb.char <- apply(comb, 2, paste, collapse = ':')
  
  # COUNTS ---------------------------------------------------------------------
  n <- tapply(values, gr, length) # group sizes
  k <- length(n) # number of groups
  N <- sum(n) # total number of observations in all groups
  
  # RANKS ----------------------------------------------------------------------
  ranks <- rank(values, ties.method = 'average') # ranks for all values
  ranks_gr <- split(ranks, gr) # ranks in each group
  r_mean <- sapply(ranks_gr, mean) # mean rank per groups
  # r_sum <- sapply(ranks_gr, sum) # sum of ranks per groups
  r_sd <- sapply(ranks_gr, sd)
  
  # CONOVER TEST STATS ---------------------------------------------------------
  R <- sapply(ranks_gr, function(i) {i^2})
  R <- sum(sapply(R, sum))
  s2 <- (1/(N-1))*(R - (N*(N+1)^2)/4)
  MS_within <- s2 * ((N-1-H) / (N-k)) # MSwithin from ANOVA calculated on ranks
  
  res <- as.data.frame(t(apply(comb, 2, conover_pair_comp, n, MS_within, r_mean)))
  colnames(res) <- c('t','SE')
  rownames(res) <- comb.char
  res$pval <- pt(abs(res$t), df = N-k, lower.tail = F) * 2 # p-value or a two-sided test
  res$p.adj <- p.adjust(res$pval, method = method)
  res$rank.cohen.d <- apply(comb, 2, cohen_rank_d, n, r_mean, r_sd)
  
  t <- t(res)[1, ]
  
  res$conover.d <- apply(comb, 2, d_Conover, n, t)
  
  # MODIFIED COHEN'S D FOR CONOVER'S RANK TESTSING -----------------------------
  return(t(res))
}


eta2_KW <- function(H, n, k) {
  # Function to calculate eta2 effect size for non-parametric Kruskal-Wallis statistics.
  # H - Kruskal-Wallis test statistics
  # n - total number of observations
  # k - number of groups
  eta2 <- (H-k-1)/(n-k)
  return(eta2)
}


d_Conover <- function(current_gr, n, t) {
  # Function to calculate Cohen's d effect size modification for non-parametric Conover post-hoc test, based on the Conover t test statistics.
  # current_groups - names of groups compared
  # t - data frame with t values from the Conover post-hoc test for a given feature
  # group_labels - group labels for each observation
  #gr <- unlist(strsplit(current_groups, split = ':', fixed = T))
  #t <- t[current_groups]
  n1<-n[current_gr][1]
  n2<-n[current_gr][2]
  t <- t[paste(as.vector(current_gr), collapse = ":")]
  # n1 <- length(which(group_labels == gr[1]))
  # n2 <- length(which(group_labels == gr[2]))
  d <- t * sqrt(1/n1 + 1/n2)
  return(as.numeric(abs(d)))
}

