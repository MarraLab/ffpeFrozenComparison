source("/projects/marralab/cayan_prj/pMedTechComparison/required.R")

args <- commandArgs(trailingOnly = TRUE)
index <- args[1]
obj_list <- readRDS(paste0("/projects/marralab/cayan_prj/pMedTechComparison/Data/batchCorrection/Objects/", index))

cfg <- list(
  cell_var    = "cellType_revised",   # cell type label column in meta.data
  batch_var   = "mode",    # FFPE vs Fresh
  ffpe_label   = "FFPE",        # exact label used for FFPE in platform column
  fresh_label  = "Fresh",       # exact label used for Fresh in platform column
  patient_var = "sample",
  n_pcs       = 30,
  dims_use    = 1:15,          # dims used for neighbors/clustering/metrics
  resolutions = seq(0.05, 1, by = 0.05) # grid to optimize NMI/ARI vs clustering
)

## ---------------- Helpers ----------------
# Safe hvgs + normalization per object
preprocess_list <- function(obj_list, assay = "RNA") {
  lapply(obj_list, function(x) {
    DefaultAssay(x) <- assay
    x <- JoinLayers(x, assay = assay)
    x <- DietSeurat(x, assay = assay, counts = TRUE, data = TRUE, scale.data = FALSE)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x
  })
}

## ----------- Integration runners -----------
run_integration <- function(obj_list, method = c("CCAIntegration","RPCAIntegration","JointPCAIntegration","HarmonyIntegration","FastMNN")) {
  method <- match.arg(method)
  obj_list <- preprocess_list(obj_list)
  
  if (method %in% c("CCAIntegration","RPCAIntegration")) {
    # Anchor-based (FindIntegrationAnchors + IntegrateData) with reduction set
    reduction_method <- switch(method,
                               CCAIntegration   = "cca",
                               RPCAIntegration  = "rpca",
                               JointPCAIntegration = "jpca")
    features <- intersect(SelectIntegrationFeatures(object.list = obj_list), Reduce(intersect, lapply(obj_list, Features)))
    obj_list <- lapply(obj_list, function(x) {
      ScaleData(x, assay = "RNA", features = features, verbose = FALSE) %>%
        RunPCA(verbose = F, features = features, layer = "scale.data", npcs = cfg$n_pcs)
    })
    anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, reduction = reduction_method, dims = cfg$dims_use)
    integrated <- IntegrateData(anchorset = anchors, k.weight = 20)
    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(integrated, verbose = FALSE)
    integrated <- RunPCA(integrated, npcs = cfg$n_pcs) # integrated PC
    integrated@reductions$emb <- integrated@reductions$pca        # unify name for scoring
    emb_name <- "emb"
    
  } else if (method == "JointPCAIntegration"){
    merged <- Reduce(merge, obj_list) %>%
      ScaleData() %>%
      RunPCA(npcs = cfg$n_pcs) %>%
      IntegrateLayers(method = JointPCAIntegration,
                      orig.reduction = "pca",
                      new.reduction  = "emb")
    emb_name <- "emb"
    integrated <- merged %>%
      JoinLayers(assay = "RNA") %>%
      NormalizeData()
  } else if (method == "HarmonyIntegration") {
    merged <- Reduce(merge, obj_list)
    merged <- ScaleData(merged, verbose = FALSE)
    merged <- RunPCA(merged, npcs = cfg$n_pcs)
    merged <- RunHarmony(merged, group.by.vars = "mode", reduction.save = "harmony", verbose = FALSE)
    emb_name <- "harmony"
    merged@reductions$emb <- merged@reductions$harmony
    integrated <- merged %>%
      JoinLayers(assay = "RNA") %>%
      NormalizeData()
    
  } else if (method == "FastMNN") {
    integrated <- RunFastMNN(object.list = obj_list)  # creates 'mnn' reduction
    emb_name <- "mnn"
    integrated@reductions$emb <- integrated@reductions$mnn
    integrated <- integrated %>%
      JoinLayers(assay = "RNA") %>%
      NormalizeData()
  }
  
  # Standard downstream graph and clustering on integrated embedding
  integrated <- FindNeighbors(integrated, reduction = emb_name, dims = 1:30)
  # choose resolution maximizing NMI vs ground-truth cell labels
  best_res <- NA; best_nmi <- -Inf; best_ari <- -Inf
  for (res in cfg$resolutions) {
    tmp <- FindClusters(integrated, resolution = res, verbose = FALSE)
    nmi <- NMI(tmp@meta.data[[cfg$cell_var]], tmp@meta.data$seurat_clusters, variant = "max")
    ari <- ARI(tmp@meta.data[[cfg$cell_var]], tmp@meta.data$seurat_clusters)
    if (nmi > best_nmi) {
      best_nmi <- nmi; best_ari <- ari; best_res <- res
      integrated <- tmp
    }
  }
  integrated@misc$best_resolution <- best_res
  integrated
}

## ----------- Metrics -----------
# Returns a 'dist' object using the fastest available backend
.get_dist <- function(X, metric = c("euclidean", "cosine", "angular", "manhattan", "hamming")) {
  metric <- match.arg(metric)

  if (requireNamespace("distances", quietly = TRUE)) {
    # 'distances' returns its own class; coerce to 'dist'
    d <- distances::distances(X, method = metric)
    return(as.dist(d))
  } else if (requireNamespace("parallelDist", quietly = TRUE)) {
    return(parallelDist::parDist(X, method = metric))
  } else {
    # stats::dist supports a subset of metrics; map unsupported ones
    method_map <- list(euclidean = "euclidean",
                       manhattan = "manhattan")
    if (!metric %in% names(method_map)) {
      warning(sprintf("Metric '%s' not supported by stats::dist; using 'euclidean'.", metric))
    }
    return(dist(X, method = method_map[[metric]] %||% "euclidean"))
  }
}

# ASW for cell type on a given embedding
score_asw_cell <- function(obj, emb = "emb",
                           metric = c("euclidean", "cosine", "angular", "manhattan", "hamming"),
                           per_type_cap = 4000,
                           max_total   = 30000) {
  metric <- match.arg(metric)

  # Extract embedding on selected dims
  emb_mat <- Embeddings(obj, emb)[, cfg$dims_use, drop = FALSE]
  labs    <- obj@meta.data[[cfg$cell_var]]
  labs_f  <- factor(labs)

  # --- Subsample per cell type, then cap total ---
  set.seed(42)
  idx <- unlist(lapply(split(seq_len(nrow(emb_mat)), labs_f), function(ix) {
    if (length(ix) <= per_type_cap) ix else sample(ix, per_type_cap)
  }))
  if (length(idx) > max_total) idx <- sample(idx, max_total)

  X    <- emb_mat[idx, , drop = FALSE]
  labS <- labs_f[idx]

  # Pairwise distances via best available backend
  d <- .get_dist(X, metric = metric)

  # Silhouette (higher is better for within–cell-type cohesion)
  sil <- cluster::silhouette(as.integer(labS), d)
  mean(sil[, "sil_width"], na.rm = TRUE)
}

# Scaled ASW for batch (platform) following scIB: compute per-cell-type abs(sil), invert & average
score_asw_batch_scaled <- function(obj, emb = "emb",
                                   metric = c("euclidean", "cosine", "angular", "manhattan", "hamming"),
                                   per_type_cap = 4000,
                                   min_cells_ct = 3) {
  metric <- match.arg(metric)

  # embedding on selected dims
  emb_mat <- Embeddings(obj, emb)[, cfg$dims_use, drop = FALSE]
  labs_ct <- factor(obj@meta.data[[cfg$cell_var]])
  batch   <- factor(obj@meta.data[[cfg$batch_var]])

  ct_levels <- levels(labs_ct)
  scores <- numeric(length(ct_levels))
  scores[] <- NA_real_

  set.seed(42)
  for (j in seq_along(ct_levels)) {
    ix <- which(labs_ct == ct_levels[j])

    # need at least min_cells_ct cells in this cell type
    if (length(ix) < min_cells_ct) next

    # cap per cell type to keep pairwise distances tractable
    if (length(ix) > per_type_cap) ix <- sample(ix, per_type_cap)

    bj <- droplevels(batch[ix])
    # need at least 2 batch levels to compute batch silhouette
    if (nlevels(bj) < 2) next

    Xj <- emb_mat[ix, , drop = FALSE]
    dj <- .get_dist(Xj, metric = metric)

    # compute silhouette safely
    sil <- tryCatch({
      cluster::silhouette(as.integer(bj), dj)
    }, error = function(e) {
      # e.g., "incorrect number of dimensions" or other edge cases
      NULL
    })

    if (is.null(sil)) next

    # some degenerate cases may not yield a 2-D matrix; check before indexing
    sil_mat <- tryCatch({
      as.matrix(sil)
    }, error = function(e) NULL)

    if (is.null(sil_mat) || !is.matrix(sil_mat) || !"sil_width" %in% colnames(sil_mat)) next

    scores[j] <- mean(1 - abs(sil_mat[, "sil_width"]), na.rm = TRUE)
  }

  # Average across cell types that produced valid scores
  mean(scores, na.rm = TRUE)
}

# NMI & ARI against ground truth cell types from current clustering
score_nmi_ari <- function(obj) {
  nmi <- NMI(obj@meta.data[[cfg$cell_var]], obj@meta.data$seurat_clusters, variant = "max")
  ari <- ARI(obj@meta.data[[cfg$cell_var]], obj@meta.data$seurat_clusters)
  c(NMI = nmi, ARI = ari)
}

# PCR: principal component regression of platform on PCs (pre vs post) and scaled reduction
pcr_variance_explained <- function(pc_scores, batch_vec, pc_var) {
  # For each PC: lm(PC ~ batch) and take R^2; weight by variance explained of PC
  r2 <- sapply(seq_len(ncol(pc_scores)), function(j) {
    # batch_vec as numeric {0,1} for 2 platforms
    y <- pc_scores[, j]
    m <- summary(lm(y ~ batch_vec))
    max(0, m$r.squared)
  })
  # weighted average by PC variance contribution
  sum(r2 * pc_var) / sum(pc_var)
}

score_pcr_comparison <- function(obj_pre, obj_post, emb_post = "emb") {
  # PRE: PCA on RNA@data
  DefaultAssay(obj_pre) <- "RNA"
  obj_pre <- FindVariableFeatures(obj_pre)
  obj_pre <- ScaleData(obj_pre, verbose = FALSE)
  obj_pre <- RunPCA(obj_pre, npcs = cfg$n_pcs, verbose = FALSE)
  pre_scores <- Embeddings(obj_pre, "pca")
  pre_var    <- obj_pre@reductions$pca@stdev^2
  
  # POST: use existing embedding; if not PCA, run PCA on integrated assay to mirror scIB flexibility
  post_scores <- Embeddings(obj_post, emb_post)
  # For non-PCA embeddings, approximate equal variance per dimension
  post_var <- rep(1, ncol(post_scores))
  
  # batch numeric
  batch <- obj_post@meta.data[[cfg$batch_var]]
  batch_num <- as.numeric(factor(batch)) - 1
  
  r2_pre  <- pcr_variance_explained(pre_scores[, cfg$dims_use, drop=FALSE],  batch_num, pre_var[cfg$dims_use])
  r2_post <- pcr_variance_explained(post_scores[, cfg$dims_use, drop=FALSE], batch_num, post_var[cfg$dims_use])
  
  # Scaled difference (0 = no change, 1 = full removal)
  score <- ifelse(r2_pre > 0, max(0, (r2_pre - r2_post) / r2_pre), 0)
  list(r2_pre = r2_pre, r2_post = r2_post, pcr_scaled = score)
}

## ----------- End-to-end: scenario × method -----------
# obj_list: list of Seurat objects (FFPE & fresh)

evaluate_all <- function(obj_list,
                         methods   = c("CCAIntegration","RPCAIntegration","JointPCAIntegration","HarmonyIntegration","FastMNN")) {
    results <- list()
    # Pre-integration merged object for PCR baseline
    pre_merged <- Reduce(merge, obj_list)
    pre_merged <- JoinLayers(pre_merged, assay = "RNA") %>%
        NormalizeData(assay = "RNA")

    fresh <- subset(pre_merged, mode == "3' scRNA")
    ffpe <- subset(pre_merged, mode == "FFPE Flex")

    for (m in methods) {
        message("  Method: ", m)

        if(m == "FastMNN"){
          obj_list2 <- obj_list
          obj_list <- list(fresh, ffpe)
        }
        integrated <- tryCatch(run_integration(obj_list, m), error = function(e) {message(e); NULL})
        if (is.null(integrated)) next
        
        emb_name <- if ("emb" %in% names(integrated@reductions)) "emb" else "pca"
        
        asw_cell  <- score_asw_cell(integrated, emb = emb_name)
        asw_batch <- score_asw_batch_scaled(integrated, emb = emb_name)   # scaled
        nm_ari    <- score_nmi_ari(integrated)
        pcr       <- score_pcr_comparison(pre_merged, integrated, emb_post = emb_name)
        
        results[[length(results)+1]] <- tibble(
            method   = m,
            ASW_cell = asw_cell,
            ASW_batch_scaled = asw_batch,
            NMI = nm_ari["NMI"],
            ARI = nm_ari["ARI"],
            PCR_R2_pre  = pcr$r2_pre,
            PCR_R2_post = pcr$r2_post,
            PCR_scaled  = pcr$pcr_scaled,
            best_resolution = integrated@misc$best_resolution
        )
    }
    out <- bind_rows(results) %>%
        mutate(freshSamples = paste(unique(fresh$sample), collapse = "|"),
               ffpeSample = paste(unique(ffpe$sample), collapse = "|"),
               freshCellType = paste(unique(fresh$cellType_revised), collapse = "|"),
               ffpeCellType = paste(unique(ffpe$cellType_revised), collapse = "|"),
               minNumCells = min(unlist(lapply(obj_list2, ncol))))
    write_csv(out, paste0("/projects/marralab/scratch/cayan/techComp/Outputs/", index, "_out.csv"))
}

evaluate_all(obj_list)

