source("required.R")
allSeu <- readRDS("Data/allCells_noBatchCorrection.RDS")

samples <- unique(allSeu$sample)
cellType <- unique(allSeu$cellType_revised)

allSeu_modeSep <- SplitObject(allSeu, split.by = "mode")

# ---- Setup ----
suppressPackageStartupMessages({
  library(future)
  library(future.apply)
})

n_workers <- max(1, min(10, parallel::detectCores() - 1))
plan(multisession, workers = n_workers)

# Optional: increase allowed size for exported globals (big Seurat objects)
options(future.globals.maxSize = 100 * 1024^3)  # 2 GB; raise if needed

# Ensure output directory exists
out_dir <- "Data/batchCorrection/Objects"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# partial overlap
with_progress({
  p <- progressor(along = 1:100)

  future_lapply(
    X = 1:100,
    FUN = function(i, max_attempts = 50, min_cells_per_sample = 2) {
      # Defensive: load packages inside worker
      suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
      })

      attempt <- 0L
      success <- FALSE
      fresh_obj <- NULL
      ffpe_obj  <- NULL

      repeat {
        attempt <- attempt + 1L

        # --- Random selections (size chosen uniformly from 1..length) ---
        ct_fresh <- sample(cellType, size = sample(seq_len(length(cellType)), size = 1))
        ct_ffpe  <- sample(cellType, size = sample(seq_len(length(cellType)), size = 1))
        samp_f   <- sample(samples,  size = sample(seq_len(length(samples)),  size = 1))
        samp_p   <- sample(samples,  size = sample(seq_len(length(samples)),  size = 1))

        # --- Require overlap in either cell types OR samples ---
        cond_ok <- (length(intersect(ct_fresh, ct_ffpe)) > 0) ||
                   (length(intersect(samp_f,  samp_p))   > 0)

        if (!cond_ok) {
          if (attempt >= max_attempts) break
          next
        }

        # --- Subset FRESH; catch "No cells found" and retry ---
        fresh_obj <- tryCatch({
          subset(
            allSeu_modeSep[[1]],
            subset = sample %in% samp_f & cellType_revised %in% ct_fresh
          )
        }, error = function(e) {
          NULL
        })
        if (is.null(fresh_obj) || ncol(fresh_obj) == 0L) {
          if (attempt >= max_attempts) break
          next
        }

        # Keep FRESH samples with >= min_cells_per_sample cells
        tabCells_fresh <- table(fresh_obj$sample) %>%
          as.data.frame() %>%
          dplyr::filter(Freq >= min_cells_per_sample)
        fresh_obj <- subset(fresh_obj, sample %in% tabCells_fresh$Var1)
        if (ncol(fresh_obj) == 0L) {
          if (attempt >= max_attempts) break
          next
        }

        # --- Subset FFPE; catch "No cells found" and retry ---
        ffpe_obj <- tryCatch({
          subset(
            allSeu_modeSep[[2]],
            subset = sample %in% samp_p & cellType_revised %in% ct_ffpe
          )
        }, error = function(e) {
          NULL
        })
        if (is.null(ffpe_obj) || ncol(ffpe_obj) == 0L) {
          if (attempt >= max_attempts) break
          next
        }

        # Keep FFPE samples with >= min_cells_per_sample cells
        tabCells_ffpe <- table(ffpe_obj$sample) %>%
          as.data.frame() %>%
          dplyr::filter(Freq >= min_cells_per_sample)
        ffpe_obj <- subset(ffpe_obj, sample %in% tabCells_ffpe$Var1)
        if (ncol(ffpe_obj) == 0L) {
          if (attempt >= max_attempts) break
          next
        }

        # --- Valid scenario achieved ---
        success <- TRUE
        break
      }

      # --- Build output (split by sample) ---
      out_list <- list()
      index    <- paste("incompleteOverlap", i, sep = "_")
      out_file <- file.path(out_dir, paste0(index, ".RDS"))

      if (success) {
        out_list <- c(
          SplitObject(fresh_obj, split.by = "sample"),
          SplitObject(ffpe_obj,  split.by = "sample")
        )
      } else {
        warning(sprintf(
          "Iter %d: failed to sample a valid scenario after %d attempts; saving empty list.",
          i, attempt
        ))
      }

      # --- Save (non-fatal on error) ---
      tryCatch({
        saveRDS(out_list, out_file)
        p(message = sprintf("Saved %s (attempts=%d, success=%s)",
                            basename(out_file), attempt, success))
      }, error = function(e) {
        warning(sprintf("Failed to save iteration %d: %s", i, conditionMessage(e)))
      })

      invisible(out_file)
    },
    future.seed = TRUE  # reproducible random sampling across workers
  )
})

# total overlap
future_lapply(
    X = 1:25,
    FUN = function(i) {
        # Load required pkgs in worker (defensive)
        suppressPackageStartupMessages({
        library(Seurat)
        })

        numSamples <- sample(seq(1, length(samples)), size = 1)
        whichSamples <- sample(samples, numSamples)
        toIntegrate <- lapply(allSeu_modeSep, function(x) subset(x, sample %in% whichSamples) %>% SplitObject(split.by = "sample"))
        index <- paste("totalOverlap", i, sep = "_")

        # Build output (split by sample, then combine the lists)
        saveRDS(unlist(toIntegrate), paste0("Data/batchCorrection/Objects/", index, ".RDS"))
        out_list <- unlist(toIntegrate)
        out_file <- file.path(out_dir, paste0(index, ".RDS"))

        # Wrap save in tryCatch to avoid hard stops
        tryCatch({
            saveRDS(out_list, out_file)
        }, error = function(e) {
            warning(sprintf("Failed to save iteration %d: %s", i, conditionMessage(e)))
        })
    },
    future.seed = TRUE  # reproducible random sampling    future.seed = TRUE  # reproducible random sampling across workers
)

# no overlap
with_progress({
  p <- progressor(along = 1:25)

  future_lapply(
    X = 1:25,
    FUN = function(i, max_attempts = 50, min_cells_per_sample = 2) {
      suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
      })

      attempt <- 0L
      success <- FALSE
      fresh   <- NULL
      ffpe    <- NULL

      repeat {
        attempt <- attempt + 1L

        ## --- Random selections (size ∈ [1..length]) ---
        cellTypes_fresh <- sample(cellType, size = sample(seq_len(length(cellType)), size = 1))
        samples_fresh   <- sample(samples,  size = sample(seq_len(length(samples)),  size = 1))

        ## --- Subset FRESH; catch "No cells found" and retry ---
        fresh <- tryCatch({
          subset(
            allSeu_modeSep[[1]],
            subset = sample %in% samples_fresh & cellType_revised %in% cellTypes_fresh
          )
        }, error = function(e) {
          # Typical message: "No cells found"; re-sample
          NULL
        })

        if (is.null(fresh) || ncol(fresh) == 0L) {
          if (attempt >= max_attempts) break
          next
        }

        # Keep samples with >= min_cells_per_sample cells
        tabCells_fresh <- table(fresh$sample) %>%
          as.data.frame() %>%
          dplyr::filter(Freq >= min_cells_per_sample)
        fresh <- subset(fresh, sample %in% tabCells_fresh$Var1)

        if (ncol(fresh) == 0L) { # after filtering it might become empty; retry
          if (attempt >= max_attempts) break
          next
        }

        ## --- FFPE "no-overlap" indices relative to FRESH selections ---
        sampleInd <- which(!(allSeu_modeSep[[2]]$sample %in% samples_fresh))
        cellTypeInd <- which(!(allSeu_modeSep[[2]]$cellType_revised %in% cellTypes_fresh))

        indices <- intersect(sampleInd,cellTypeInd)
        cells <- colnames(allSeu_modeSep[[2]])[indices]

        if (length(cells) < 2L) {  # need at least 2 cells to proceed
          if (attempt >= max_attempts) break
          next
        }

        ## --- Subset FFPE safely; catch "No cells found" and retry ---
        ffpe <- tryCatch({
          subset(allSeu_modeSep[[2]], cells = cells)
        }, error = function(e) {
          NULL
        })

        if (is.null(ffpe) || ncol(ffpe) == 0L) {
          if (attempt >= max_attempts) break
          next
        }

        # Keep FFPE samples with >= min_cells_per_sample cells
        tabCells_ffpe <- table(ffpe$sample) %>%
          as.data.frame() %>%
          dplyr::filter(Freq >= min_cells_per_sample)
        ffpe <- subset(ffpe, sample %in% tabCells_ffpe$Var1)

        # Validate both sides have cells after filtering
        if (ncol(ffpe) > 0L && ncol(fresh) > 0L) {
          success <- TRUE
          break
        }

        if (attempt >= max_attempts) break
      }

      ## --- Assemble and save result ---
      out_list <- list()
      index    <- paste("noOverlap", i, sep = "_")
      out_file <- file.path(out_dir, paste0(index, ".RDS"))

      if (success) {
        out_list <- c(
          SplitObject(fresh, split.by = "sample"),
          SplitObject(ffpe,  split.by = "sample")
        )
      } else {
        warning(sprintf(
          "Iter %d: failed to build a valid scenario after %d attempts; saving empty list.",
          i, attempt
        ))
      }

      tryCatch({
        saveRDS(out_list, out_file)
        p(message = sprintf("Saved %s (attempts=%d, success=%s)",
                            basename(out_file), attempt, success))
      }, error = function(e) {
        warning(sprintf("Failed to save iteration %d: %s", i, conditionMessage(e)))
      })

      invisible(out_file)
    },
    future.seed = TRUE  # reproducible random sampling across workers
  )
})