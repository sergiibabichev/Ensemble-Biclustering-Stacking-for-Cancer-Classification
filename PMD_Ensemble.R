# ============================================================
# PMD (PMA) ensemble with MV/time/coverage logging + consensus
# (fixed PMD call, robust u/v -> biclusters, sanity test)
# ============================================================

suppressPackageStartupMessages({
  library(PMA)
})

# ---------- common utils ----------
load_full_matrix <- function(file_path) {
  df <- read.csv(file_path, check.names = FALSE)
  expr <- t(as.matrix(df[, -ncol(df)]))
  rownames(expr) <- gsub("^X", "", rownames(expr))
  storage.mode(expr) <- "double"
  if (anyNA(expr) || any(!is.finite(expr))) stop("NA/Inf in matrix.")
  expr
}

mv_for_block <- function(X, rows, cols) {
  if (length(rows) < 2 || length(cols) < 2) return(NA_real_)
  B <- X[rows, cols, drop = FALSE]
  if (!all(is.finite(B))) return(NA_real_)
  rmean <- rowMeans(B); cmean <- colMeans(B); gmean <- mean(B)
  R <- B - rmean %*% matrix(1,1,ncol(B)) - matrix(1,nrow(B),1) %*% t(cmean) + gmean
  tot_var <- stats::var(as.numeric(B)); if (!is.finite(tot_var) || tot_var==0) return(NA_real_)
  res_var <- stats::var(as.numeric(R))
  as.numeric(res_var / tot_var)
}

mean_mv_for_structure <- function(X, biclusters, min_genes=10, min_samples=10) {
  if (length(biclusters)==0) return(list(mean_mv=NA_real_, mv_per_cluster=numeric(0)))
  mv_vals <- numeric(0)
  for (bc in biclusters) {
    rows <- bc$rows; cols <- bc$cols
    if (length(rows) > min_genes && length(cols) > min_samples) {
      mv <- tryCatch(mv_for_block(X, rows, cols), error=function(e) NA_real_)
      if (is.finite(mv)) mv_vals <- c(mv_vals, mv)
    }
  }
  if (!length(mv_vals)) return(list(mean_mv=NA_real_, mv_per_cluster=numeric(0)))
  list(mean_mv=mean(mv_vals), mv_per_cluster=mv_vals)
}

# ---------- data ----------
X_full <- load_full_matrix("cancer_Combine_filtered.csv")
G <- nrow(X_full); S <- ncol(X_full)
cat(sprintf("PMD: data %d genes x %d samples\n", G, S))

# ---------- ensemble params ----------
B_runs <- 60
g_min <- 6000; g_max <- 7000
set.seed(2025)

# PMD optimal params (Table 4.2)
sumabsu_base <- 40
sumabsv_base <- 5
K_pmd   <- 10  # number of bicluster components per run

min_genes_bc <- 10
min_samples_bc <- 10

H_hits <- matrix(0L, nrow=G, ncol=S)
D_gene <- integer(G)
run_log <- data.frame(run=integer(0), n_genes=integer(0), n_biclusters=integer(0),
                      mv_mean=numeric(0), runtime_sec=numeric(0), coverage_frac=numeric(0))
mv_per_cluster_dump <- vector("list", B_runs)

# ---------- helpers: build biclusters from u/v ----------
# participation-ratio based top-mass selector (коли підтримка порожня/занадто широка)
select_support_by_effsize <- function(w) {
  a <- abs(w); s <- sum(a)
  if (!is.finite(s) || s == 0) return(integer(0))
  p <- a / s
  n_eff <- max(1, round(1 / sum(p^2)))  # participation ratio
  ord <- order(p, decreasing = TRUE)
  ord[seq_len(min(n_eff, length(ord)))]
}

# PMD -> список бікластерів: якщо u/v мають нулі — беремо ненульові;
# якщо вся колона майже щільна або пуста — fallback на select_support_by_effsize
pmd_to_bics <- function(U, V, min_g=10, min_s=10) {
  K <- min(ncol(U), ncol(V)); if (K == 0) return(list())
  out <- list(); nb <- 0L
  for (k in seq_len(K)) {
    u <- U[,k]; v <- V[,k]
    r <- which(abs(u) > 0); c <- which(abs(v) > 0)
    if (length(r) == 0 || length(r) == nrow(U)) r <- select_support_by_effsize(u)
    if (length(c) == 0 || length(c) == ncol(V)) c <- select_support_by_effsize(v)
    if (length(r) >= min_g && length(c) >= min_s) {
      nb <- nb + 1L; out[[nb]] <- list(rows = sort(r), cols = sort(c))
    }
  }
  out
}

# ---------- one run ----------
run_pmd_once <- function() {
  n_g <- sample(g_min:g_max, 1)
  genes_idx <- sample.int(G, n_g, replace=FALSE)
  Xg <- X_full[genes_idx, , drop=FALSE]
  
  # PMA::PMD має обмеження L1: sumabsu <= sqrt(#rows), sumabsv <= sqrt(#cols)
  maxU <- sqrt(nrow(Xg)); maxV <- sqrt(ncol(Xg))
  su <- min(sumabsu_base, maxU)
  sv <- min(sumabsv_base, maxV)
  
  t0 <- proc.time()[3]
  fit <- tryCatch(
    PMA::PMD(
      x       = as.matrix(Xg),
      type    = "standard",   # ФІКСУЄМО режим
      K       = K_pmd,
      sumabsu = su,
      sumabsv = sv,
      center  = TRUE,         # центровані колонки (працює для експресії)
      trace   = FALSE
    ),
    error = function(e) { message("PMD error: ", e$message); NULL }
  )
  t1 <- proc.time()[3]
  runtime <- as.numeric(t1 - t0)
  
  if (is.null(fit) || is.null(fit$u) || is.null(fit$v)) {
    return(list(ok=FALSE, genes_idx=genes_idx, n_biclusters=0, mv_mean=NA_real_,
                mv_per_cluster=numeric(0), runtime=runtime,
                mask_hits_local=integer(0), coverage=0))
  }
  
  # U: genes x K, V: samples x K
  U <- as.matrix(fit$u)
  V <- as.matrix(fit$v)
  
  bics <- pmd_to_bics(U, V, min_g = min_genes_bc, min_s = min_samples_bc)
  if (!length(bics)) {
    return(list(ok=FALSE, genes_idx=genes_idx, n_biclusters=0, mv_mean=NA_real_,
                mv_per_cluster=numeric(0), runtime=runtime,
                mask_hits_local=integer(0), coverage=0))
  }
  
  mv_res <- mean_mv_for_structure(Xg, bics, min_genes_bc, min_samples_bc)
  mv_mean <- mv_res$mean_mv
  
  mask_local <- matrix(0L, nrow=nrow(Xg), ncol=ncol(Xg))
  for (k in seq_along(bics)) {
    r <- bics[[k]]$rows; c <- bics[[k]]$cols
    if (length(r) && length(c)) mask_local[r,c] <- mask_local[r,c] + 1L
  }
  coverage <- mean(mask_local > 0L)
  
  list(ok=TRUE, genes_idx=genes_idx, n_biclusters=length(bics),
       mv_mean=mv_mean, mv_per_cluster=mv_res$mv_per_cluster,
       runtime=runtime, mask_hits_local=mask_local, coverage=coverage)
}

# --- Quick sanity test on a random 3000-gene subset ---
set.seed(2025)
test_idx <- sample.int(G, min(3000L, G), replace = FALSE)
X_test <- X_full[test_idx, , drop = FALSE]
maxU_test <- sqrt(nrow(X_test)); maxV_test <- sqrt(ncol(X_test))
su_test <- min(sumabsu_base, maxU_test); sv_test <- min(sumabsv_base, maxV_test)

test_fit <- tryCatch(
  PMA::PMD(as.matrix(X_test), type="standard", K=K_pmd,
           sumabsu=su_test, sumabsv=sv_test, center=TRUE, trace=FALSE),
  error = function(e) { message("Test PMD error: ", e$message); NULL }
)
if (!is.null(test_fit) && !is.null(test_fit$u) && !is.null(test_fit$v)) {
  tb <- pmd_to_bics(as.matrix(test_fit$u), as.matrix(test_fit$v),
                    min_g=min_genes_bc, min_s=min_samples_bc)
  cov_tb <- if (length(tb)) {
    m <- matrix(0L, nrow=nrow(X_test), ncol=ncol(X_test))
    for (i in seq_along(tb)) m[tb[[i]]$rows, tb[[i]]$cols] <- 1L
    mean(m > 0L)
  } else 0
  cat(sprintf("Sanity test biclusters: %d, cov=%.3f\n", length(tb), cov_tb))
} else {
  cat("Sanity test: PMD fit is NULL\n")
}

# ---------- ensemble loop ----------
cat(sprintf("Starting PMD ensemble: %d runs; gene subset [%d,%d]\n", B_runs, g_min, g_max))
for (b in seq_len(B_runs)) {
  cat(sprintf("PMD run %d/%d ... ", b, B_runs))
  out <- run_pmd_once()
  D_gene[out$genes_idx] <- D_gene[out$genes_idx] + 1L
  if (out$ok) {
    H_hits[out$genes_idx, ] <- H_hits[out$genes_idx, ] + (out$mask_hits_local > 0L)
    mv_per_cluster_dump[[b]] <- out$mv_per_cluster
    cat(sprintf("clusters=%d, MV_mean=%.4f, time=%.1fs, cov=%.3f\n",
                out$n_biclusters, out$mv_mean, out$runtime, out$coverage))
  } else {
    mv_per_cluster_dump[[b]] <- numeric(0)
    cat(sprintf("no clusters, time=%.1fs\n", out$runtime))
  }
  run_log <- rbind(run_log, data.frame(run=b, n_genes=length(out$genes_idx),
                                       n_biclusters=out$n_biclusters, mv_mean=out$mv_mean,
                                       runtime_sec=out$runtime, coverage_frac=out$coverage))
}

# ---------- consensus ----------
D_mat <- matrix(pmax(D_gene,1L), nrow=G, ncol=S)
Consensus <- H_hits / D_mat

dir.create("ensemble_pmd", showWarnings=FALSE)
write.csv(run_log, "ensemble_pmd/PMD_ensemble_runlog.csv", row.names=FALSE)
saveRDS(mv_per_cluster_dump, "ensemble_pmd/PMD_mv_per_cluster_list.rds")
saveRDS(Consensus, "ensemble_pmd/PMD_consensus_matrix.rds")

cat("\nPMD summary:\n")
cat(sprintf("Valid runs: %d/%d\n", sum(is.finite(run_log$mv_mean)), B_runs))
cat(sprintf("Mean MV_mean: %.4f\n", mean(run_log$mv_mean, na.rm=TRUE)))
cat(sprintf("Mean runtime (s): %.1f\n", mean(run_log$runtime_sec, na.rm=TRUE)))
cat(sprintf("Mean coverage: %.3f\n", mean(run_log$coverage_frac, na.rm=TRUE)))
