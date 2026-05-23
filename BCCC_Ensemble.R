# ============================================================
# Ensembling BCCC over bootstrap-like gene subsamples
# Logs: MV_mean per run, runtime, cluster counts
# Consensus: frequency of membership per cell, normalized by gene sampling count
# ============================================================

# ----------------------------
# 0) Libs
# ----------------------------
suppressPackageStartupMessages({
  library(biclust)
  library(Matrix)   # на випадок, якщо захочете перейти на sparse
})

# ----------------------------
# 1) Завантаження повної матриці (genes x samples)
#    (ваш формат: samples в стовпцях CSV, останній стовпець — ярлик; ми його відкидаємо)
# ----------------------------
load_full_matrix <- function(file_path) {
  df <- read.csv(file_path, check.names = FALSE)
  expr <- t(as.matrix(df[, -ncol(df)]))
  rownames(expr) <- gsub("^X", "", rownames(expr))
  storage.mode(expr) <- "double"
  if (anyNA(expr) || any(!is.finite(expr))) {
    stop("❌ У матриці є NA/Inf. Приберіть/імпутуйте перед запуском.")
  }
  expr
}

X_full <- load_full_matrix("cancer_Combine_filtered.csv")  # genes x samples
G <- nrow(X_full); S <- ncol(X_full)
cat(sprintf("Data loaded: %d genes x %d samples\n", G, S))

# ----------------------------
# 2) Ядро метрики (MV)
# ----------------------------
mv_for_block <- function(X, rows, cols) {
  if (length(rows) < 2 || length(cols) < 2) return(NA_real_)
  B <- X[rows, cols, drop = FALSE]
  if (!all(is.finite(B))) return(NA_real_)
  rmean <- rowMeans(B); cmean <- colMeans(B); gmean <- mean(B)
  # подвійне центрювання:
  # R = B - rmean - cmean + gmean
  R <- B - rmean %*% matrix(1, 1, ncol(B)) - matrix(1, nrow(B), 1) %*% t(cmean) + gmean
  tot_var <- stats::var(as.numeric(B))
  if (!is.finite(tot_var) || tot_var == 0) return(NA_real_)
  res_var <- stats::var(as.numeric(R))
  as.numeric(res_var / tot_var)
}

mean_mv_for_structure <- function(X, biclusters, min_genes = 10, min_samples = 10) {
  if (length(biclusters) == 0) return(list(mean_mv = NA_real_, mv_per_cluster = numeric(0)))
  mv_vals <- numeric(0)
  for (bc in biclusters) {
    rows <- bc$rows; cols <- bc$cols
    if (length(rows) > min_genes && length(cols) > min_samples) {
      mv <- tryCatch(mv_for_block(X, rows, cols), error = function(e) NA_real_)
      if (is.finite(mv)) mv_vals <- c(mv_vals, mv)
    }
  }
  if (length(mv_vals) == 0) return(list(mean_mv = NA_real_, mv_per_cluster = numeric(0)))
  list(mean_mv = mean(mv_vals), mv_per_cluster = mv_vals)
}

biclust_to_list <- function(bc_res) {
  out <- vector("list", bc_res@Number)
  for (k in seq_len(bc_res@Number)) {
    rows <- which(bc_res@RowxNumber[, k])
    cols <- which(bc_res@NumberxCol[k, ])
    out[[k]] <- list(rows = rows, cols = cols)
  }
  out
}

# ----------------------------
# 3) Параметри ансамблювання
# ----------------------------
B_runs     <- 60                # кількість прогонів
g_min      <- 6000              # мінімум генів у прогоні
g_max      <- 7000              # максимум генів у прогоні
set.seed(2025)                  # відтворюваність

# оптимальні параметри BCCC згідно табл. 4.2
bccc_delta <- 1.2
bccc_alpha <- 1.0
bccc_number <- 50               # максимальна кількість бікластерів (налаштовується)

# мінімальні розміри бікластерів для MV
min_genes_bc   <- 10
min_samples_bc <- 10

# ----------------------------
# 4) Акумулюючі структури
# ----------------------------
# H: скільки разів клітинка (gene i, sample j) потрапила у будь-який бікластер (сума лічильників)
# D: скільки разів ген i був включений у підвибірку (для нормування по рядках)
H_hits <- matrix(0L, nrow = G, ncol = S)   # integer матриця (економніше за numeric)
D_gene <- integer(G)                       # лічильник включень гена

run_log <- data.frame(
  run = integer(0),
  n_genes = integer(0),
  n_biclusters = integer(0),
  mv_mean = numeric(0),
  runtime_sec = numeric(0),
  coverage_frac = numeric(0),    # частка клітин у підматриці, що потрапили до ≥1 бікластеру
  stringsAsFactors = FALSE
)

mv_per_cluster_dump <- vector("list", B_runs)  # для збереження розподілу MV по кластерах

# ----------------------------
# 5) Один прогін
# ----------------------------
run_bccc_once <- function(run_id) {
  n_g <- sample(g_min:g_max, size = 1)
  genes_idx <- sample.int(G, size = n_g, replace = FALSE)
  Xg <- X_full[genes_idx, , drop = FALSE]
  
  t0 <- proc.time()[[3]]
  bc_res <- tryCatch(
    biclust(Xg, method = BCCC(), delta = bccc_delta, alpha = bccc_alpha, number = bccc_number),
    error = function(e) NULL
  )
  t1 <- proc.time()[[3]]
  runtime <- as.numeric(t1 - t0)
  
  if (is.null(bc_res) || bc_res@Number == 0) {
    return(list(
      ok = FALSE,
      genes_idx = genes_idx,
      n_biclusters = 0,
      mv_mean = NA_real_,
      mv_per_cluster = numeric(0),
      runtime = runtime,
      mask_hits_local = integer(0),
      coverage = 0
    ))
  }
  
  # список бікластерів
  bics <- biclust_to_list(bc_res)
  
  # MV
  mv_res <- mean_mv_for_structure(Xg, bics, min_genes = min_genes_bc, min_samples = min_samples_bc)
  mv_mean <- mv_res$mean_mv
  
  # побудова маски підпрогона (підматриця)
  mask_local <- matrix(0L, nrow = nrow(Xg), ncol = ncol(Xg))
  for (k in seq_len(bc_res@Number)) {
    r <- which(bc_res@RowxNumber[, k])
    c <- which(bc_res@NumberxCol[k, ])
    if (length(r) > 0 && length(c) > 0) {
      mask_local[r, c] <- mask_local[r, c] + 1L
    }
  }
  coverage <- mean(mask_local > 0)
  
  list(
    ok = TRUE,
    genes_idx = genes_idx,
    n_biclusters = bc_res@Number,
    mv_mean = mv_mean,
    mv_per_cluster = mv_res$mv_per_cluster,
    runtime = runtime,
    mask_hits_local = mask_local,
    coverage = coverage
  )
}

# ----------------------------
# 6) Основний цикл ансамблю
#    (за бажанням можна паралелити через parallel::mclapply)
# ----------------------------
cat(sprintf("Starting ensemble: %d runs; gene subset in [%d, %d]\n", B_runs, g_min, g_max))

for (b in seq_len(B_runs)) {
  cat(sprintf(">>> Run %d/%d ... ", b, B_runs))
  out <- run_bccc_once(b)
  
  # оновлюємо лічильники
  D_gene[out$genes_idx] <- D_gene[out$genes_idx] + 1L
  
  if (out$ok) {
    # “піднімаємо” локальну маску на глобальні індекси
    H_hits[out$genes_idx, ] <- H_hits[out$genes_idx, ] + (out$mask_hits_local > 0L)
    mv_per_cluster_dump[[b]] <- out$mv_per_cluster
    cat(sprintf("clusters=%d, MV_mean=%.4f, time=%.1fs, coverage=%.3f\n",
                out$n_biclusters, out$mv_mean, out$runtime, out$coverage))
  } else {
    mv_per_cluster_dump[[b]] <- numeric(0)
    cat(sprintf("no clusters, time=%.1fs\n", out$runtime))
  }
  
  # лог рядок
  run_log <- rbind(run_log, data.frame(
    run = b,
    n_genes = length(out$genes_idx),
    n_biclusters = out$n_biclusters,
    mv_mean = out$mv_mean,
    runtime_sec = out$runtime,
    coverage_frac = out$coverage,
    stringsAsFactors = FALSE
  ))
}

# ----------------------------
# 7) Консенсус: нормалізація по кількості включень гена
#    F[i, j] = H_hits[i, j] / D_gene[i]
# ----------------------------
D_mat <- matrix(pmax(D_gene, 1L), nrow = G, ncol = S)  # запобігаємо /0
Consensus <- H_hits / D_mat

# ----------------------------
# 8) Збереження результатів
# ----------------------------
dir.create("ensemble_bccc", showWarnings = FALSE)

write.csv(run_log, file = "ensemble_bccc/BCCC_ensemble_runlog.csv", row.names = FALSE)
saveRDS(mv_per_cluster_dump, file = "ensemble_bccc/BCCC_mv_per_cluster_list.rds")
saveRDS(Consensus, file = "ensemble_bccc/BCCC_consensus_matrix.rds")

# коротке резюме
cat("\n=== Ensemble summary ===\n")
cat(sprintf("Valid runs: %d of %d\n", sum(is.finite(run_log$mv_mean)), B_runs))
cat(sprintf("Mean MV_mean over valid runs: %.4f\n", mean(run_log$mv_mean, na.rm = TRUE)))
cat(sprintf("Median MV_mean: %.4f\n", median(run_log$mv_mean, na.rm = TRUE)))
cat(sprintf("Mean runtime (s): %.1f\n", mean(run_log$runtime_sec, na.rm = TRUE)))
cat(sprintf("Mean coverage: %.3f\n", mean(run_log$coverage_frac, na.rm = TRUE)))
cat("Consensus matrix saved to ensemble_bccc/BCCC_consensus_matrix.rds\n")
