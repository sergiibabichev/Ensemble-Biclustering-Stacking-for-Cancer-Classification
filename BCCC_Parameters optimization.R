
# ----------------------------
# 1) Завантаження та відбір генів
# ----------------------------
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

# Перевірка коректності матриці
if (any(is.na(full_matrix)) || any(!is.finite(full_matrix))) {
  stop("❌ У матриці виявлені NA або Inf значення")
}

# ---------------------------
# 2) Пакети
# ---------------------------
library(biclust)
library(ggplot2)
library(viridis)  # можна прибрати, якщо використовуєте тільки scale_fill_viridis_c() з ggplot2

# ---------------------------
# 3) Ядро метрики (універсальні функції)
# ---------------------------

# MV для однієї підматриці
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

# ---------------------------
# 4) Сітка параметрів
# ---------------------------
params_grid <- expand.grid(
  delta = seq(1.2, 2.0, by = 0.2),
  alpha = seq(1.0, 1.4, by = 0.2)
)

# ---------------------------
# 5) Grid Search для BCCC
# ---------------------------
results <- data.frame()

for (i in 1:nrow(params_grid)) {
  delta <- params_grid$delta[i]
  alpha <- params_grid$alpha[i]
  cat(sprintf(">>> [%d/%d] delta = %.2f, alpha = %.2f\n", i, nrow(params_grid), delta, alpha))
  
  res_bccc <- tryCatch({
    biclust(full_matrix, method = BCCC(),
            delta = delta,
            alpha = alpha,
            number = 50)
  }, error = function(e) NULL)
  
  if (!is.null(res_bccc)) {
    bics <- biclust_to_list(res_bccc)
    mv_res <- mean_mv_for_structure(full_matrix, bics, min_genes = 10, min_samples = 10)
    bccc_mv <- mv_res$mean_mv
    results <- rbind(results, data.frame(delta = delta, alpha = alpha, MV_BCCC = bccc_mv))
    cat(sprintf("  → MV = %.4f\n", bccc_mv))
  } else {
    results <- rbind(results, data.frame(delta = delta, alpha = alpha, MV_BCCC = NA_real_))
    cat("  ⚠️ Помилка при побудові моделі\n")
  }
}

# ---------------------------
# 6) Найкраща комбінація + перевірка, що є валідні MV
# ---------------------------
if (any(!is.na(results$MV_BCCC))) {
  best_row <- results[which.min(results$MV_BCCC), ]
  cat(sprintf("\n✅ Найкраща комбінація: delta = %.2f, alpha = %.2f → MV = %.4f\n",
              best_row$delta, best_row$alpha, best_row$MV_BCCC))
} else {
  cat("\n⚠️ Усі MV = NA. Перевірте фільтри min_genes/min_samples або параметри BCCC.\n")
}

write.csv(results, "BCCC_results.csv", row.names = FALSE)


# ---------------------------
# Функція для побудови масок
# ---------------------------
build_mask_matrix <- function(bc_result, nrow, ncol) {
  mask <- matrix(0, nrow = nrow, ncol = ncol)
  for (k in 1:bc_result@Number) {
    rows <- bc_result@RowxNumber[, k]
    cols <- bc_result@NumberxCol[k, ]
    mask[rows, cols] <- mask[rows, cols] + 1
  }
  return(mask)
}











