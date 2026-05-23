# =========================
# 0) Завантаження даних
# =========================
load_data <- function(file_path, top_n) {
  df <- read.csv(file_path, check.names = FALSE)
  expr <- t(as.matrix(df[, -ncol(df)]))
  gene_vars <- apply(expr, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:top_n]
  mat <- expr[top_genes, ]
  rownames(mat) <- gsub("^X", "", rownames(mat))
  mat
}

full_matrix <- load_data("cancer_Combine_filtered.csv", top_n = 6000)

# Перевірка на NA/Inf
if (anyNA(full_matrix) || any(!is.finite(full_matrix))) {
  stop("❌ У матриці виявлені NA або Inf значення")
}


# ---------------------------
# 3) Метрика якості (універсальна MV)
# ---------------------------

# MV для однієї підматриці (подвійне центрювання через sweep)
mv_for_block <- function(X, rows, cols) {
  if (length(rows) < 2 || length(cols) < 2) return(NA_real_)
  B <- X[rows, cols, drop = FALSE]
  if (!all(is.finite(B))) return(NA_real_)
  rmean <- rowMeans(B); cmean <- colMeans(B); gmean <- mean(B)
  R <- B - rmean %*% matrix(1, 1, ncol(B)) - matrix(1, nrow(B), 1) %*% t(cmean) + gmean
  tot_var <- stats::var(as.numeric(B))
  if (!is.finite(tot_var) || tot_var == 0) return(NA_real_)
  res_var <- stats::var(as.numeric(R))
  as.numeric(res_var / tot_var)
}

# Середнє MV по структурі (список бікластерів)
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

# Конвертер Biclust -> список індексів
biclust_to_list <- function(bc_res) {
  out <- vector("list", bc_res@Number)
  for (k in seq_len(bc_res@Number)) {
    rows <- which(bc_res@RowxNumber[, k])
    cols <- which(bc_res@NumberxCol[k, ])
    out[[k]] <- list(rows = rows, cols = cols)
  }
  out
}

# ===== PMA::PMD (Sparse SVD) =====
# install.packages("PMA")  # якщо потрібно
# ===== PMA::PMD (Sparse SVD) — grid search без ручної бінаризації =====
# install.packages("PMA")  # якщо потрібно
library(PMA)

# --- ВАША метрика (використовую як є; заберіть, якщо вже визначена вище) ---
# MV для однієї підматриці

# --- Автовибір підтримки для випадків без розрідження (без ручних порогів) ---
select_support_by_effsize <- function(w) {
  a <- abs(w); s <- sum(a)
  if (!is.finite(s) || s == 0) return(integer(0))
  p <- a / s
  n_eff <- max(1, round(1 / sum(p^2)))  # participation ratio
  ord <- order(p, decreasing = TRUE)
  ord[seq_len(min(n_eff, length(ord)))]
}

# PMD -> список бікластерів: ненульова підтримка; якщо її нема/вся — auto top-N
pmd_to_list_nz <- function(pmd) {
  U <- pmd$u; V <- pmd$v
  if (is.null(U) || is.null(V)) return(list())
  K <- ncol(U)
  out <- vector("list", K)
  for (k in seq_len(K)) {
    rows <- which(U[, k] != 0)
    cols <- which(V[, k] != 0)
    if (length(rows) == 0 || length(rows) == nrow(U)) rows <- select_support_by_effsize(U[, k])
    if (length(cols) == 0 || length(cols) == nrow(V)) cols <- select_support_by_effsize(V[, k])
    out[[k]] <- list(rows = sort(rows), cols = sort(cols))
  }
  out
}

# ---- Параметри гріду ----
maxU <- sqrt(nrow(full_matrix))
maxV <- sqrt(ncol(full_matrix))

params_pmd <- expand.grid(
  sumabsu = c(5, 10, 20, 40),
  sumabsv = c(5, 10, 20, 40)
)

K_PMD <- 10
results_pmd <- data.frame()

# ---- Перебір ----
set.seed(123)
for (i in seq_len(nrow(params_pmd))) {
  su <- params_pmd$sumabsu[i]; sv <- params_pmd$sumabsv[i]
  cat(sprintf("PMD   [%03d/%03d] sumabsu=%.1f sumabsv=%.1f  (K=%d)\n",
              i, nrow(params_pmd), su, sv, K_PMD))
  
  mv_mean <- tryCatch({
    pmd <- PMA::PMD(
      as.matrix(full_matrix),
      type    = "standard",   # важливо зафіксувати саме стандартний режим
      K       = K_PMD,
      sumabsu = min(su, maxU),
      sumabsv = min(sv, maxV),
      center  = FALSE,
      trace   = FALSE
    )
    if (is.null(pmd$u) || is.null(pmd$v)) stop("PMD returned NULL u/v")
    
    bics <- pmd_to_list_nz(pmd)  # ← без порогів/квантилів
    mean_mv_for_structure(full_matrix, bics, min_genes = 10, min_samples = 10)$mean_mv
  }, error = function(e) { cat("  ⚠️", e$message, "\n"); NA_real_ })
  
  results_pmd <- rbind(results_pmd, data.frame(sumabsu = su, sumabsv = sv, K = K_PMD, MV = mv_mean))
  if (is.finite(mv_mean)) cat(sprintf("  → MV = %.4f\n", mv_mean))
}

# ---- Підсумок ----
if (any(!is.na(results_pmd$MV))) {
  best_pmd <- results_pmd[which.min(results_pmd$MV), ]
  cat(sprintf("✅ PMD best: sumabsu=%.1f  sumabsv=%.1f  K=%d  MV=%.4f\n",
              best_pmd$sumabsu, best_pmd$sumabsv, best_pmd$K, best_pmd$MV))
} else {
  cat("⚠️ PMD: усі MV=NA — збільшіть sumabsu/sumabsv, зменшіть вимоги min_genes/min_samples або перевірте, що u/v не порожні.\n")
}

write.csv(results_pmd, "PMD_results.csv", row.names = FALSE)

