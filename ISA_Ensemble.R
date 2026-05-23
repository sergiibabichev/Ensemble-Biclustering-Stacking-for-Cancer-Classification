# =========================================================
# ISA ensemble on FULL matrix (no prefilter), with MV, time,
# per-run coverage, and consensus accumulation over full universe
# + core→expansion to boost coverage, optional soft retry
# =========================================================

suppressPackageStartupMessages({
  library(isa2)
})

# --- 1) Завантаження повної матриці (гени x зразки) ---
load_full_matrix <- function(file_path) {
  df <- read.csv(file_path, check.names = FALSE)
  expr <- t(as.matrix(df[, -ncol(df)]))
  rownames(expr) <- gsub("^X", "", rownames(expr))
  storage.mode(expr) <- "double"
  if (anyNA(expr) || any(!is.finite(expr))) stop("NA/Inf in matrix.")
  expr
}

full_matrix <- load_full_matrix("cancer_Combine_filtered.csv")
stopifnot(is.matrix(full_matrix))

# --- 2) Метрики MV ---
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

# --- 3) Запуск ISA (коректні аргументи для isa2::isa) ---
isa_run_once <- function(Xg, thr_row, thr_col, no_seeds) {
  fit <- isa2::isa(
    Xg,
    thr.row   = thr_row,
    thr.col   = thr_col,
    no.seeds  = no_seeds,
    direction = c("updown","updown")
  )
  if (is.null(fit$rows) || is.null(fit$columns)) return(list(biclusters = list(), k = 0))
  K <- min(ncol(fit$rows), ncol(fit$columns))
  if (is.null(K) || K == 0) return(list(biclusters = list(), k = 0))
  bics <- vector("list", K)
  for (k in seq_len(K)) {
    rows <- which(fit$rows[, k] != 0)
    cols <- which(fit$columns[, k] != 0)
    bics[[k]] <- list(rows = rows, cols = cols)
  }
  list(biclusters = bics, k = K)
}

# --- 4) Параметри ISA та ансамблю ---
set.seed(2025)

isa_params <- list(
  thr_row   = 0.30,
  thr_col   = 3.00,
  no_seeds  = 100
)

# опційні “м’які” пороги, якщо локальне cov зовсім мале
isa_soft <- list(
  thr_row = 0.25,
  thr_col = 2.75
)

ensemble_cfg <- list(
  n_boot         = 60,     # кількість бутстрепів
  genes_per_boot = 6500,   # ~6–7k генів як обговорювали
  min_genes      = 10,
  min_samples    = 10
)

# --- Expansion params (post-ISA growth) ---
EXPAND_ROWS   <- TRUE
EXPAND_COLS   <- TRUE
TAU_ROW       <- 0.25   # corr threshold gene ↔ профіль зразків модуля
TAU_COL       <- 0.25   # corr threshold sample ↔ профіль генів модуля
MAX_ADD_ROWS  <- 1500
MAX_ADD_COLS  <- 500
MV_TOLERANCE  <- 0.10   # допустиме відносне погіршення MV (+10%); 0 = заборонити
COV_RETRY_THR <- 0.02   # якщо локальне покриття нижче — пробуємо м'які пороги

# --- 5) Акумулятори ---
G_all <- nrow(full_matrix)
S_all <- ncol(full_matrix)

consensus_mask <- matrix(0L, nrow = G_all, ncol = S_all)

boot_summary <- data.frame(
  run = integer(0),
  n_genes = integer(0),
  k_modules = integer(0),
  mv_mean = numeric(0),
  coverage_frac = numeric(0),
  elapsed_sec = numeric(0)
)

out_dir <- "isa_ensemble_full_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 5.1) Допоміжні для expansion ---
safe_cor <- function(a, b) {
  a <- a - mean(a); b <- b - mean(b)
  sa <- sqrt(sum(a*a)); sb <- sqrt(sum(b*b))
  if (sa == 0 || sb == 0) return(0)
  sum(a*b)/(sa*sb)
}

expand_bicluster <- function(X, rows, cols,
                             tau_row = 0.25, tau_col = 0.25,
                             max_add_rows = 1500, max_add_cols = 500,
                             mv_tolerance = 0.10,
                             min_genes = 10, min_samples = 10,
                             expand_rows = TRUE, expand_cols = TRUE) {
  if (length(rows) < min_genes || length(cols) < min_samples)
    return(list(rows = rows, cols = cols))
  
  sub <- X[rows, cols, drop = FALSE]
  rcent <- colMeans(sub)  # профіль по зразках
  ccent <- rowMeans(sub)  # профіль по генах
  
  add_rows <- integer(0)
  if (expand_rows) {
    cand_r <- setdiff(seq_len(nrow(X)), rows)
    cors_r <- vapply(cand_r, function(i) safe_cor(X[i, cols], rcent), numeric(1))
    ord_r  <- order(cors_r, decreasing = TRUE)
    keep_r <- cand_r[ord_r][cors_r[ord_r] >= tau_row]
    if (length(keep_r) > max_add_rows) keep_r <- keep_r[seq_len(max_add_rows)]
    add_rows <- keep_r
  }
  
  add_cols <- integer(0)
  if (expand_cols) {
    cand_c <- setdiff(seq_len(ncol(X)), cols)
    cors_c <- vapply(cand_c, function(j) safe_cor(X[rows, j], ccent), numeric(1))
    ord_c  <- order(cors_c, decreasing = TRUE)
    keep_c <- cand_c[ord_c][cors_c[ord_c] >= tau_col]
    if (length(keep_c) > max_add_cols) keep_c <- keep_c[seq_len(max_add_cols)]
    add_cols <- keep_c
  }
  
  new_rows <- sort(unique(c(rows, add_rows)))
  new_cols <- sort(unique(c(cols, add_cols)))
  
  if (mv_tolerance < Inf) {
    mv_core <- mv_for_block(X, rows, cols)
    mv_new  <- mv_for_block(X, new_rows, new_cols)
    if (is.finite(mv_core) && is.finite(mv_new) && mv_new > (1 + mv_tolerance) * mv_core) {
      tmp_rows <- sort(unique(c(rows, add_rows))); tmp_cols <- cols
      mv_try1  <- mv_for_block(X, tmp_rows, tmp_cols)
      tmp_rows2 <- rows; tmp_cols2 <- sort(unique(c(cols, add_cols)))
      mv_try2   <- mv_for_block(X, tmp_rows2, tmp_cols2)
      
      choices <- list(
        list(r=rows,      c=cols,      mv=mv_core),
        list(r=tmp_rows,  c=tmp_cols,  mv=mv_try1),
        list(r=tmp_rows2, c=tmp_cols2, mv=mv_try2),
        list(r=new_rows,  c=new_cols,  mv=mv_new)
      )
      valid <- vapply(choices, function(z) is.finite(z$mv) && z$mv <= (1+mv_tolerance)*mv_core, logical(1))
      if (any(valid)) {
        idx <- which(valid)[which.min(vapply(choices[valid], `[[`, numeric(1), "mv"))]
        return(list(rows=choices[[idx]]$r, cols=choices[[idx]]$c))
      } else {
        return(list(rows=rows, cols=cols))
      }
    }
  }
  list(rows=new_rows, cols=new_cols)
}

# --- 6) Ансамблевий цикл ---
for (b in seq_len(ensemble_cfg$n_boot)) {
  cat(sprintf("ISA bootstrap %d/%d ...\n", b, ensemble_cfg$n_boot))
  t0 <- proc.time()
  
  # 6.1) Випадкова підмножина генів з ПОВНОГО універсуму
  genes_subset_idx <- sample.int(G_all, size = ensemble_cfg$genes_per_boot, replace = FALSE)
  Xg <- full_matrix[genes_subset_idx, , drop = FALSE]
  
  # 6.2) Запуск ISA (базові пороги)
  isa_res <- tryCatch(
    isa_run_once(Xg, isa_params$thr_row, isa_params$thr_col, isa_params$no_seeds),
    error = function(e) { message("ISA error: ", e$message); NULL }
  )
  if (is.null(isa_res)) {
    elapsed <- (proc.time() - t0)[["elapsed"]]
    boot_summary <- rbind(
      boot_summary,
      data.frame(run = b, n_genes = nrow(Xg), k_modules = 0,
                 mv_mean = NA_real_, coverage_frac = NA_real_, elapsed_sec = elapsed)
    )
    next
  }
  
  # 6.3) Дорощування кожного бікластеру (core → expansion)
  bics_expanded <- isa_res$biclusters
  if (length(bics_expanded) > 0) {
    for (k in seq_along(bics_expanded)) {
      bc <- bics_expanded[[k]]
      if (length(bc$rows) && length(bc$cols)) {
        bics_expanded[[k]] <- expand_bicluster(
          Xg, bc$rows, bc$cols,
          tau_row = TAU_ROW, tau_col = TAU_COL,
          max_add_rows = MAX_ADD_ROWS, max_add_cols = MAX_ADD_COLS,
          mv_tolerance = MV_TOLERANCE,
          min_genes = ensemble_cfg$min_genes,
          min_samples = ensemble_cfg$min_samples,
          expand_rows = EXPAND_ROWS, expand_cols = EXPAND_COLS
        )
      }
    }
  }
  
  # 6.4) Локальна маска + локальне покриття (після expansion)
  local_mask <- matrix(0L, nrow = nrow(Xg), ncol = ncol(Xg))
  if (length(bics_expanded) > 0) {
    for (k in seq_along(bics_expanded)) {
      r <- bics_expanded[[k]]$rows
      c <- bics_expanded[[k]]$cols
      if (length(r) && length(c)) local_mask[r, c] <- local_mask[r, c] + 1L
    }
  }
  coverage <- mean(local_mask > 0L)
  
  # 6.4b) Якщо cov дуже мале — однократний "м'який" повтор із нижчими порогами
  if (coverage < COV_RETRY_THR) {
    isa_res2 <- tryCatch(
      isa_run_once(Xg, thr_row = isa_soft$thr_row, thr_col = isa_soft$thr_col, no_seeds = isa_params$no_seeds),
      error = function(e) NULL
    )
    if (!is.null(isa_res2) && length(isa_res2$biclusters)) {
      b2 <- lapply(isa_res2$biclusters, function(bc)
        expand_bicluster(Xg, bc$rows, bc$cols,
                         tau_row=TAU_ROW, tau_col=TAU_COL,
                         max_add_rows=MAX_ADD_ROWS, max_add_cols=MAX_ADD_COLS,
                         mv_tolerance=MV_TOLERANCE,
                         min_genes=ensemble_cfg$min_genes, min_samples=ensemble_cfg$min_samples,
                         expand_rows=EXPAND_ROWS, expand_cols=EXPAND_COLS))
      local2 <- matrix(0L, nrow=nrow(Xg), ncol=ncol(Xg))
      for (k in seq_along(b2)) {
        r <- b2[[k]]$rows; c <- b2[[k]]$cols
        if (length(r) && length(c)) local2[r,c] <- local2[r,c] + 1L
      }
      # беремо максимум по клітинці (збільшує покриття)
      local_mask <- pmax(local_mask, local2)
      coverage <- mean(local_mask > 0L)
      # також оновимо список модулів для оцінки MV
      bics_expanded <- c(bics_expanded, b2)
    }
  }
  
  # 6.5) Якість за MV (на локальній підматриці) — на expanded/soft
  mv_res <- mean_mv_for_structure(
    Xg, bics_expanded,
    min_genes = ensemble_cfg$min_genes,
    min_samples = ensemble_cfg$min_samples
  )
  mv_mean <- mv_res$mean_mv
  
  # 6.6) Проєкція локальної маски у ПОВНИЙ універсум
  consensus_mask[genes_subset_idx, ] <- consensus_mask[genes_subset_idx, ] + local_mask
  
  # 6.7) Логи часу/якості
  elapsed <- (proc.time() - t0)[["elapsed"]]
  boot_summary <- rbind(
    boot_summary,
    data.frame(run = b, n_genes = nrow(Xg), k_modules = length(bics_expanded),
               mv_mean = mv_mean, coverage_frac = coverage, elapsed_sec = elapsed)
  )
  
  # збереження артефактів прогону
  saveRDS(
    list(
      genes_subset_idx = genes_subset_idx,
      biclusters = bics_expanded,
      mv_per_cluster = mv_res$mv_per_cluster,
      mv_mean = mv_mean,
      coverage_frac = coverage,
      elapsed_sec = elapsed,
      params = list(isa = isa_params, isa_soft = isa_soft,
                    expand = list(TAU_ROW=TAU_ROW, TAU_COL=TAU_COL,
                                  MAX_ADD_ROWS=MAX_ADD_ROWS, MAX_ADD_COLS=MAX_ADD_COLS,
                                  MV_TOLERANCE=MV_TOLERANCE))
    ),
    file = file.path(out_dir, sprintf("isa_run_%02d.rds", b))
  )
  
  cat(sprintf("  → modules: %d, MV_mean = %.4f, cov = %.3f, time = %.2fs\n",
              length(bics_expanded), mv_mean, coverage, elapsed))
}

# --- 7) Збереження підсумків ансамблю ---
write.csv(boot_summary, file.path(out_dir, "isa_boot_summary.csv"), row.names = FALSE)
saveRDS(consensus_mask, file.path(out_dir, "isa_consensus_mask_full.rds"))

# --- 8) Глобальні метрики покриття по консенсусу ---
global_cell_coverage <- mean(consensus_mask > 0L)
gene_coverage   <- rowMeans(consensus_mask > 0L)  # частка зразків, де ген «накритий»
sample_coverage <- colMeans(consensus_mask > 0L)  # частка генів, «накритих» у зразку

write.csv(data.frame(metric = "global_cell_coverage", value = global_cell_coverage),
          file.path(out_dir, "isa_consensus_global_coverage.csv"), row.names = FALSE)
write.csv(data.frame(gene = seq_len(nrow(consensus_mask)), coverage = gene_coverage),
          file.path(out_dir, "isa_gene_coverage.csv"), row.names = FALSE)
write.csv(data.frame(sample = seq_len(ncol(consensus_mask)), coverage = sample_coverage),
          file.path(out_dir, "isa_sample_coverage.csv"), row.names = FALSE)

cat("\n✅ ISA ensemble (FULL matrix) finished.\n")
cat(sprintf("Avg MV_mean (valid runs): %.4f\n",
            mean(boot_summary$mv_mean[is.finite(boot_summary$mv_mean)])))
cat(sprintf("Avg coverage (valid runs): %.3f\n",
            mean(boot_summary$coverage_frac[is.finite(boot_summary$coverage_frac)])))
cat(sprintf("Global consensus coverage: %.3f\n", global_cell_coverage))
cat(sprintf("Avg time per run: %.2fs\n",
            mean(boot_summary$elapsed_sec)))
