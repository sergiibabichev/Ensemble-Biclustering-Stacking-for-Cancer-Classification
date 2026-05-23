suppressPackageStartupMessages({
  library(matrixStats)
  library(data.table)
  library(stats); library(utils)
  library(biclust)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(DOSE)
  library(ggplot2)
})

## --------- PATHS ---------
DATA_CSV <- "cancer_Combine_filtered.csv"   
CONS_ISA  <- "ISA_consensus_matrix.rds"
CONS_PMD  <- "PMD_consensus_matrix.rds"
OUT_DIR   <- "stacking_output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


## --------- 1) LOADERS & UTILS ---------
load_expression <- function(file_path) {
  df <- read.csv(file_path, check.names = FALSE)
  # припущення: останній стовпчик — labels (як у ваших скриптах)
  expr <- t(as.matrix(df[, -ncol(df)]))
  rownames(expr) <- gsub("^X", "", rownames(expr))
  storage.mode(expr) <- "double"
  if (anyNA(expr) || any(!is.finite(expr))) stop("NA/Inf in expression matrix.")
  labels <- df[[ncol(df)]]
  list(X = expr, labels = labels, df_raw = df)
}

normalize_matrix01 <- function(M) {
  mx <- max(M, na.rm=TRUE)
  if (mx > 0 && is.finite(mx)) M / mx else M
}

align_masks_by_names <- function(masks) {
  # Залишаємо лише спільні імена рядків/стовпців та впорядковуємо
  rn_common <- Reduce(intersect, lapply(masks, rownames))
  cn_common <- Reduce(intersect, lapply(masks, colnames))
  if (length(rn_common) < 2 || length(cn_common) < 2)
    stop("Too few common row/col names across consensus matrices.")
  lapply(masks, function(M) M[rn_common, cn_common, drop=FALSE])
}

stack_masks <- function(masks, weights=NULL) {
  stopifnot(length(masks) >= 1)
  if (is.null(weights)) weights <- rep(1/length(masks), length(masks))
  stopifnot(length(weights) == length(masks))
  S <- 0
  for (i in seq_along(masks)) S <- S + weights[i] * masks[[i]]
  S
}

## Greedy extractor on continuous S (no binarization)
extract_biclusters <- function(S, tau=0.7, max_k=150, min_rows=10, min_cols=12, X_ref=NULL) {
  p <- nrow(S); n <- ncol(S)
  B <- list(); M <- S; iter <- 0
  while (iter < max_k && is.finite(max(M))) {
    if (max(M) < tau) break
    iter <- iter + 1
    idx <- which(M == max(M), arr.ind=TRUE)[1,]
    rows <- idx[1]; cols <- idx[2]
    changed <- TRUE
    while (changed) {
      changed <- FALSE
      if (!length(cols)) break
      rscore <- rowMeans(M[, cols, drop=FALSE])
      new_r  <- which(rscore >= tau)
      if (!all(new_r %in% rows)) { rows <- new_r; changed <- TRUE }
      if (!length(rows)) break
      cscore <- colMeans(M[rows, , drop=FALSE])
      new_c  <- which(cscore >= tau)
      if (!all(new_c %in% cols)) { cols <- new_c; changed <- TRUE }
    }
    if (length(rows) >= min_rows && length(cols) >= min_cols) {
      B[[length(B)+1]] <- list(rows=rows, cols=cols, src="STACK")
      M[rows, cols] <- -Inf
    } else {
      M[rows, cols] <- -Inf
    }
  }
  if (!is.null(X_ref) && length(B)) {
    ord <- order(sapply(B, function(b) {
      sub <- X_ref[b$rows, b$cols, drop=FALSE]
      ri <- rowMeans(sub); cj <- colMeans(sub); mu <- mean(sub)
      R  <- sub - outer(ri, rep(1, ncol(sub))) - outer(rep(1, nrow(sub)), cj) + mu
      mean(R^2)
    }))
    B <- B[ord]
  }
  B
}

## Metrics on X
mv_bicluster_X <- function(X, rows, cols) {
  sub <- X[rows, cols, drop=FALSE]
  ri <- rowMeans(sub); cj <- colMeans(sub); mu <- mean(sub)
  R  <- sub - outer(ri, rep(1, ncol(sub))) - outer(rep(1, nrow(sub)), cj) + mu
  mean(R^2)
}
msr_cheng_church_X <- mv_bicluster_X  # однакова формула для середнього квадрата залишків

safe_mean_abs_corr <- function(M, by_rows=TRUE) {
  if (by_rows) {
    if (nrow(M) < 2) return(NA_real_)
    C <- suppressWarnings(cor(t(M), use="pairwise.complete.obs"))
  } else {
    if (ncol(M) < 2) return(NA_real_)
    C <- suppressWarnings(cor(M, use="pairwise.complete.obs"))
  }
  if (nrow(C) < 2) return(NA_real_)
  mean(abs(C[upper.tri(C)]), na.rm=TRUE)
}

coverage_cells <- function(p, n, biclusters) {
  M <- matrix(FALSE, p, n)
  for (b in biclusters) if (length(b$rows) && length(b$cols)) M[b$rows, b$cols] <- TRUE
  mean(M)
}
pairwise_overlap_jaccard <- function(p, n, biclusters) {
  k <- length(biclusters); if (k < 2) return(0)
  J <- function(b1, b2){
    M1 <- matrix(FALSE, p, n); M1[b1$rows, b1$cols] <- TRUE
    M2 <- matrix(FALSE, p, n); M2[b2$rows, b2$cols] <- TRUE
    u <- sum(M1 | M2); if (u==0) return(NA_real_) else sum(M1 & M2)/u
  }
  vals <- c()
  for (i in 1:(k-1)) for (j in (i+1):k) vals <- c(vals, J(biclusters[[i]], biclusters[[j]]))
  mean(vals, na.rm=TRUE)
}

evaluate_on_X <- function(X, S, B_list, lambda=1.0, mu=0.5, target_genes=15000) {
  p <- nrow(X); n <- ncol(X); k <- length(B_list)
  MV_vals  <- vapply(B_list, function(b) mv_bicluster_X(X, b$rows, b$cols), numeric(1))
  MSR_vals <- vapply(B_list, function(b) msr_cheng_church_X(X, b$rows, b$cols), numeric(1))
  ROW_CORR <- vapply(B_list, function(b) safe_mean_abs_corr(X[b$rows, b$cols, drop=FALSE], by_rows=TRUE),  numeric(1))
  COL_CORR <- vapply(B_list, function(b) safe_mean_abs_corr(X[b$rows, b$cols, drop=FALSE], by_rows=FALSE), numeric(1))
  Coverage <- coverage_cells(p, n, B_list)
  Overlap  <- pairwise_overlap_jaccard(p, n, B_list)
  dens_S   <- vapply(B_list, function(b) if (length(b$rows) && length(b$cols)) mean(S[b$rows, b$cols]) else NA_real_, numeric(1))
  union_genes <- if (k) length(unique(unlist(lapply(B_list, `[[`, "rows")))) else 0
  
  MVmean <- if (k) mean(MV_vals, na.rm=TRUE) else NA_real_
  Obj <- MVmean + lambda * abs(union_genes - target_genes) / p + mu * (1 - Coverage)
  
  list(
    summary = data.frame(
      k = k,
      MVmean = MVmean,
      coverage_cells = Coverage,
      mean_overlap_jaccard = Overlap,
      mean_row_coherence = if (k) mean(ROW_CORR, na.rm=TRUE) else NA_real_,
      mean_col_coherence = if (k) mean(COL_CORR, na.rm=TRUE) else NA_real_,
      mean_MSR = if (k) mean(MSR_vals, na.rm=TRUE) else NA_real_,
      union_genes = union_genes,
      Obj = Obj
    ),
    per_bicluster = data.frame(
      bicluster = seq_len(k),
      n_genes   = vapply(B_list, function(b) length(b$rows), integer(1)),
      n_samples = vapply(B_list, function(b) length(b$cols), integer(1)),
      MV        = MV_vals,
      MSR       = MSR_vals,
      row_coherence = ROW_CORR,
      col_coherence = COL_CORR,
      density_on_S  = dens_S
    )
  )
}


## 1) Якщо в масок немає dimnames, але розміри збігаються з X — призначаємо імена з X
assign_dimnames_from_X <- function(M, X, label){
  if (is.null(rownames(M)) && is.null(colnames(M))) {
    stopifnot(all(dim(M) == dim(X)))
    rownames(M) <- rownames(X)
    colnames(M) <- colnames(X)
  } else {
    # якщо імена частково були — можна ще перевірити збіг із X
    if (!is.null(rownames(M))) stop(sprintf("%s unexpectedly already has rownames.", label))
    if (!is.null(colnames(M))) stop(sprintf("%s unexpectedly already has colnames.", label))
  }
  M
}


## --------- 2) LOAD DATA & CONSENSUS ---------
dat <- load_expression(DATA_CSV)
X <- dat$X
G <- nrow(X); Smp <- ncol(X)

## 0) Переконаймось, що X уже завантажений і має імена
stopifnot(nrow(X) == 18564, ncol(X) == 6344)

CM_ISA  <- readRDS(CONS_ISA)
CM_PMD  <- readRDS(CONS_PMD)

# Переводимо всі консенсуси у [0,1]
CM_ISA  <- normalize_matrix01(CM_ISA)
CM_PMD  <- normalize_matrix01(CM_PMD)

CM_ISA  <- assign_dimnames_from_X(CM_ISA,  X, "ISA")
CM_PMD  <- assign_dimnames_from_X(CM_PMD,  X, "PMD")

# Вирівнюємо маски за спільними іменами
masks <- list(ISA=CM_ISA, PMD=CM_PMD)

weights <- rep(1/length(masks), length(masks))  # рівні ваги
S <- stack_masks(unname(masks), weights=weights)
cat(sprintf("S range: [%.4f, %.4f]\n", min(S), max(S)))

## «загострення» S: посилюємо стабільні клітини
gamma <- 1.5
S2 <- S^gamma

## --------- 4) HYPERPARAM TUNING (tau by quantiles) + GREEDY EXTRACTION ---------
min_genes   <- 10
min_samples <- 12
max_k       <- 150
target_genes <- 12000
lambda <- 2    # штраф за відхилення від цільового U
mu     <- 0.6    # штраф за низьке покриття

tau_quantiles <- c(0.8, 0.85, 0.9, 0.95)
tau_grid <- as.numeric(quantile(S2, probs = tau_quantiles, na.rm=TRUE))
tau_grid <- sort(unique(tau_grid))  # на випадок збігів

tuning_results <- list()
for (tau in tau_grid) {
  B_greedy <- extract_biclusters(S2, tau=tau, max_k=max_k,
                                 min_rows=min_genes, min_cols=min_samples, X_ref=S2)
  evalX <- evaluate_on_X(X, S2, B_greedy, lambda=lambda, mu=mu, target_genes=target_genes)
  tuning_results[[length(tuning_results)+1]] <- list(
    tau = tau, q = NA_real_,
    summary = evalX$summary, per = evalX$per_bicluster, B = B_greedy
  )
  cat(sprintf("[tau=%.6f] k=%d | MVmean=%.5g | cov=%.3f | U=%d | Obj=%.4f\n",
              tau,
              evalX$summary$k,
              evalX$summary$MVmean,
              evalX$summary$coverage_cells,
              evalX$summary$union_genes,
              evalX$summary$Obj))
}

# вибір найкращого за мінімумом Obj
objs <- sapply(tuning_results, function(z) z$summary$Obj)
best_i <- which.min(objs)
best <- tuning_results[[best_i]]
B_stack <- best$B
summary_X <- best$summary
per_X <- best$per

cat("\n=== BEST (by Obj) ===\n"); print(summary_X)
cat(sprintf("Chosen tau = %.6f\n", best$tau))

# збереження таблиць
f_sum <- file.path(OUT_DIR, "STACK_short_best_summary_on_X.csv")
f_per <- file.path(OUT_DIR, "STACK_short_best_biclusters_on_X.csv")
write.csv(summary_X, f_sum, row.names=FALSE)
write.csv(per_X,   f_per, row.names=FALSE)
cat(sprintf("Saved: %s\nSaved: %s\n", f_sum, f_per))

## --------- 5) UNION GENES & ID MAPPING (to ENTREZ) ---------
union_gene_idx <- if (length(B_stack)) sort(unique(unlist(lapply(B_stack, `[[`, "rows")))) else integer(0)
union_gene_ids <- rownames(X)[union_gene_idx]

## Union gene IDs вже у форматі ENTREZID
union_entrez <- as.character(union_gene_ids)
union_entrez <- union_entrez[!is.na(union_entrez) & union_entrez != ""]
union_entrez <- unique(union_entrez)

cat(sprintf("Union genes (unique rows): %d (ENTREZ usable: %d)\n",
            length(union_gene_ids), length(union_entrez)))


## =========================
## CONSISTENT ENRICHMENT for STACK (ENTREZID, Homo sapiens)
## =========================

# =========================
# Enrichment (Homo sapiens, ENTREZID): KEGG + GO (BP/CC/MF), p/q = 0.05
# =========================
do_enrichment_hsa <- function(genes_entrez,
                              bg_entrez   = NULL,        # optional background
                              p_cutoff    = 0.05,
                              q_cutoff    = 0.05,
                              min_gs_size = 10,          # мін. розмір генного набору
                              max_gs_size = 5000,        # макс. розмір генного набору
                              out_dir     = OUT_DIR) {
  
  # підготовка вхідних генів (ENTREZID)
  genes_entrez <- unique(as.character(na.omit(genes_entrez)))
  if (!length(genes_entrez)) { message("Empty gene list; skipping."); return(NULL) }
  if (!is.null(bg_entrez))    bg_entrez <- unique(as.character(na.omit(bg_entrez)))
  
  # ---------- GO: три онтології ----------
  run_go <- function(ont) {
    tryCatch(
      enrichGO(gene          = genes_entrez,
               OrgDb         = org.Hs.eg.db,
               keyType       = "ENTREZID",
               ont           = ont,                 # "BP" / "CC" / "MF"
               pAdjustMethod = "BH",
               pvalueCutoff  = p_cutoff,
               qvalueCutoff  = q_cutoff,
               universe      = bg_entrez,
               minGSSize     = min_gs_size,
               maxGSSize     = max_gs_size,
               readable      = TRUE),
      error = function(e) NULL
    )
  }
  
  ego_bp <- run_go("BP")
  ego_cc <- run_go("CC")
  ego_mf <- run_go("MF")
  
  # об'єднаємо GO в одну табличку для зручності
  go_bind <- function(x, label) {
    if (is.null(x)) return(NULL)
    df <- as.data.frame(x)
    if (!nrow(df)) return(NULL)
    df$ONTOLOGY <- label
    df
  }
  go_all_df <- do.call(rbind, list(
    go_bind(ego_bp, "BP"),
    go_bind(ego_cc, "CC"),
    go_bind(ego_mf, "MF")
  ))
  
  if (!is.null(go_all_df) && nrow(go_all_df)) {
    write.csv(go_all_df, file.path(out_dir, "STACK_enrichment_GO_all_BP_CC_MF.csv"),
              row.names = FALSE)
  }
  
  # ---------- KEGG (людина) ----------
  # clusterProfiler очікує ENTREZID у вигляді строки, organism = "hsa"
  ekk <- tryCatch(
    enrichKEGG(gene          = genes_entrez,
               organism      = "hsa",
               pAdjustMethod = "BH",
               pvalueCutoff  = p_cutoff,
               qvalueCutoff  = q_cutoff,
               universe      = bg_entrez,
               minGSSize     = min_gs_size,
    ),
    error = function(e) NULL
  )
  if (!is.null(ekk)) {
    dfk <- as.data.frame(ekk)
    if (nrow(dfk)) write.csv(dfk, file.path(out_dir, "enrichment_KEGG.csv"),
                             row.names = FALSE)
  }
  
  list(GO = list(BP = ego_bp, CC = ego_cc, MF = ego_mf), KEGG = ekk)
}


# === ВИКЛИК ===

bg_ids <- tryCatch(unique(as.character(rownames(X))), error = function(e) NULL)

enr <- tryCatch(
  do_enrichment_hsa(union_gene_ids, bg_entrez = bg_ids,
                    p_cutoff = 0.05, q_cutoff = 0.05),
  error = function(e) { message(e$message); NULL }
)


# =========================
# KEGG enrichment + export table + dotplot
# =========================
run_kegg_with_plots <- function(genes_entrez,
                                bg_entrez   = NULL,
                                p_cutoff    = 0.05,
                                q_cutoff    = 0.05,
                                top_n       = 20,
                                out_dir     = OUT_DIR) {
  
  
  # Фіксуємо організм як Homo sapiens
  genes_entrez <- unique(as.character(na.omit(genes_entrez)))
  if (!length(genes_entrez)) stop("Gene list is empty.")
  
  if (!is.null(bg_entrez)) bg_entrez <- unique(as.character(na.omit(bg_entrez)))
  
  ekk <- enrichKEGG(gene          = genes_entrez,
                    organism      = "hsa",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = p_cutoff,
                    qvalueCutoff  = q_cutoff,
                    universe      = bg_entrez)
  
  if (is.null(ekk) || nrow(as.data.frame(ekk)) == 0) {
    warning("No enriched KEGG terms found.")
    return(NULL)
  }
  
  # Повний KEGG результат
  df_all <- as.data.frame(ekk)
  write.csv(df_all, file.path(out_dir, "STACK_short_KEGG_enrichment_full.csv"), row.names = FALSE)
  
  # Top-N терміни за p.adjust
  df_top <- df_all[order(df_all$p.adjust), ][1:min(top_n, nrow(df_all)), ]
  write.csv(df_top, file.path(out_dir, sprintf("STACK_short_KEGG_enrichment_top%d.csv", top_n)),
            row.names = FALSE)
  
  # Dotplot (clusterProfiler)
  suppressPackageStartupMessages(library(enrichplot))
  p <- dotplot(ekk, showCategory = top_n) +
    ggplot2::ggtitle(sprintf("STACK: KEGG Pathway Enrichment (Top %d)", top_n)) +
    ggplot2::theme_minimal()
  
  # Збереження графіків
  ggplot2::ggsave(file.path(out_dir, sprintf("STACK_short_KEGG_dotplot_top%d.png", top_n)),
                  plot = p, width = 8, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(out_dir, sprintf("STACK_shortKEGG_dotplot_top%d.pdf", top_n)),
                  plot = p, width = 8, height = 6)
  
  return(list(full = df_all, top = df_top, plot = p))
}

# === Виклик ===
bg_ids <- tryCatch(unique(as.character(rownames(X))), error = function(e) NULL)

kegg_res <- run_kegg_with_plots(union_gene_ids,
                                bg_entrez = bg_ids,
                                p_cutoff = 0.05,
                                q_cutoff = 0.05,
                                top_n = 20)

df <- read.csv(DATA_CSV, check.names = FALSE)
labels <- as.vector(df[, ncol(df)])
df_subset <- df[, union_gene_ids]
df_transform <- cbind(df_subset, Class = labels)

write.csv(df_transform, "STACK_short_subset_unique_genes.csv", row.names = FALSE)























