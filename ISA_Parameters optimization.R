# ===== ISA2 =====

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

X <- as.matrix(full_matrix)             # ваша початкова матриця
storage.mode(X) <- "double"
X[!is.finite(X)] <- NA

## 1) Орієнтація: гени — рядки, зразки — стовпці
# якщо у вас навпаки — розкоментуйте:
# X <- t(X)

## 2) Лог і двоетапна стандартизація
if (all(X >= 0, na.rm=TRUE) && max(X, na.rm=TRUE) > 50) {
  X <- log1p(X)
}
# z-score по рядках, потім по стовпцях
X <- t(scale(t(X), center=TRUE, scale=TRUE))
X <- scale(X, center=TRUE, scale=TRUE)
X[!is.finite(X)] <- 0

## 3) Приберемо проблемні рядки/стовпці (занадто багато NA або нульова дисперсія)
keep_r <- apply(X, 1, function(v) sd(v, na.rm=TRUE) > 0)
keep_c <- apply(X, 2, function(v) sd(v, na.rm=TRUE) > 0)
row_map <- which(keep_r); col_map <- which(keep_c)
Xn <- X[keep_r, keep_c, drop=FALSE]



# ---------------------------
# 3) Метрика якості (універсальна MV)
# ---------------------------

# MV для однієї підматриці (подвійне центрювання через sweep)
mv_for_block <- function(X, rows, cols) {
  if (length(rows) < 2 || length(cols) < 2) return(NA_real_)
  B <- X[rows, cols, drop = FALSE]
  if (!all(is.finite(B))) return(NA_real_)
  R <- sweep(sweep(B, 1, rowMeans(B), "-"), 2, colMeans(B), "-") + mean(B)
  tot_var <- stats::var(as.numeric(B))
  if (!is.finite(tot_var) || tot_var == 0) return(NA_real_)
  res_var <- stats::var(as.numeric(R))
  as.numeric(res_var / tot_var)
}

# Середній MV по структурі (список бікластерів)
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



# install.packages("isa2")  # якщо потрібно
library(isa2)

# Адаптер: ISA → список бікластерів
isa_to_list_map <- function(fit, row_map, col_map, min_genes=10, min_samples=10){
  out <- list()
  if (is.null(fit$rows) || is.null(fit$columns)) return(out)
  K <- if (is.matrix(fit$rows)) ncol(fit$rows) else length(fit$rows)
  for (k in seq_len(K)) {
    r <- which(fit$rows[,k] != 0)
    c <- which(fit$columns[,k] != 0)
    if (length(r) >= min_genes && length(c) >= min_samples) {
      out[[length(out)+1]] <- list(rows = row_map[r], cols = col_map[c])
    }
  }
  out
}

count_modules <- function(fit, min_genes=5, min_samples=5){
  if (is.null(fit$rows) || is.null(fit$columns)) return(0)
  K <- if (is.matrix(fit$rows)) ncol(fit$rows) else length(fit$rows)
  n <- 0
  for (k in seq_len(K)) {
    if (sum(fit$rows[,k]!=0) >= min_genes && sum(fit$columns[,k]!=0) >= min_samples) n <- n+1
  }
  n
}

# Сітка порогів (налаштуйте за потреби)
thr_row_grid <- c(seq(0.3, 0.9, by=0.2), seq(1.0, 3.0, by=0.25))
thr_col_grid <- c(seq(0.3, 0.9, by=0.2), seq(1.0, 3.0, by=0.25))
seed_grid    <- c(100, 200, 400, 800)   # адаптивно збільшуємо

results_isa <- data.frame(thr_row=numeric(), thr_col=numeric(), seeds=integer(),
                          n_modules=integer(), MV=as.numeric())



## 6) Пошук: спершу знайдемо хоч якісь модулі (послабимо фільтр), потім порахуємо MV
## 6) Пошук: спершу знайдемо хоч якісь модулі (послабимо фільтр 5x5), потім MV на 10x10
for (tr in thr_row_grid) for (tc in thr_col_grid) {
  cat(sprintf("ISA try thr.row=%.2f thr.col=%.2f\n", tr, tc))
  
  used_seeds <- 100L
  fit <- isa(Xn, thr.row = tr, thr.col = tc, no.seeds = used_seeds,
             direction = c("updown","updown"))
  n_mod <- count_modules(fit, 5, 5)
  cat(sprintf("  seeds=%3d -> modules>=5x5: %d\n", used_seeds, n_mod))
  
  if (n_mod == 0) {
    for (ns in c(200L, 400L, 800L)) {
      fit2 <- isa(Xn, thr.row = tr, thr.col = tc, no.seeds = ns,
                  direction = c("updown","updown"))
      n_mod2 <- count_modules(fit2, 5, 5)
      cat(sprintf("  seeds=%3d -> modules>=5x5: %d\n", ns, n_mod2))
      if (n_mod2 > 0) { 
        fit <- fit2; n_mod <- n_mod2; used_seeds <- ns
        break
      }
    }
  }
  
  bics <- isa_to_list_map(fit, row_map, col_map, min_genes = 10, min_samples = 10)
  mv   <- if (length(bics)) mean_mv_for_structure(full_matrix, bics, 10, 10)$mean_mv else NA_real_
  
  results_isa <- rbind(results_isa, data.frame(
    thr_row = tr, thr_col = tc, seeds = used_seeds, n_modules = n_mod, MV = mv
  ))
}


## 7) Підсумок
ok <- which(is.finite(results_isa$MV))
if (length(ok)) {
  best <- results_isa[ok[which.min(results_isa$MV[ok])], ]
  cat(sprintf("✅ ISA best: thr.row=%.2f thr.col=%.2f seeds=%d n_mod=%d MV=%.4f\n",
              best$thr_row, best$thr_col, best$seeds, best$n_modules, best$MV))
} else {
  cat("⚠️ ISA: досі нема валідних MV. Спробуйте зменшити min_genes/min_samples до 5/5 або розширити thr до [0.2, 3.5].\n")
}