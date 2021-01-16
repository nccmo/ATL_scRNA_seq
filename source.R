# Most of these functions are modified from functions in satijalab / seurat (GitHub; https://github.com/satijalab/seurat)
# License of Seurat: GNU General Public License v3.0

# subset v3 barcode
library(matrixStats)
library(future)
library(pbapply)

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package.installed], collapse = ', '),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}



WilcoxDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  j <- seq_len(length.out = length(x = cells.1))
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  overflow.check <- ifelse(
    test = is.na(x = suppressWarnings(length(x = data.use[1, ]) * length(x = data.use[1, ]))),
    yes = FALSE,
    no = TRUE
  )
  limma.check <- PackageCheck("limma", error = FALSE)
  if (limma.check[1] && overflow.check) {
    p_val <- my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, statistics = data.use[x, ])), 1))
        #return(min(2 * min(limma::rankSumTestWithCorrelation(index = seq_len(rowSums(!is.na(data.use[,cells.1]))[[x]]),
        #                                                    statistics = data.use[x,!is.na(data.use[x,])])),
        #          1))
      }
    )
  } else {
    if (getOption('Seurat.limma.wilcox.msg', TRUE) && overflow.check) {
      message(
        "For a more efficient implementation of the Wilcoxon Rank Sum Test,",
        "\n(default method for FindMarkers) please install the limma package",
        "\n--------------------------------------------",
        "\ninstall.packages('BiocManager')",
        "\nBiocManager::install('limma')",
        "\n--------------------------------------------",
        "\nAfter installation of limma, Seurat will automatically use the more ",
        "\nefficient implementation (no further action necessary).",
        "\nThis message will be shown once per session"
      )
      options(Seurat.limma.wilcox.msg = FALSE)
    }
    group.info <- data.frame(row.names = c(cells.1, cells.2))
    group.info[cells.1, "group"] <- "Group1"
    group.info[cells.2, "group"] <- "Group2"
    group.info[, "group"] <- factor(x = group.info[, "group"])
    data.use <- data.use[, rownames(x = group.info), drop = FALSE]
    p_val <- my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        return(wilcox.test(data.use[x, ] ~ group.info[, "group"], ...)$p.value)
      }
    )
  }
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}

FindMarkers_na_removal <- function(
  object,
  slot = "data",
  counts = numeric(),
  cells.1 = NULL, 
  cells.2 = NULL, 
  ident.1 = NULL,
  ident.2 = NULL,
  features = NULL,
  reduction = NULL,
  logfc.threshold = 0.15,
  test.use = 'wilcox',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  ...
) {
  DefaultAssay(object) <- "ADT"
  if (is.null(features)){
    features <- rownames(x = object)
  } else{
    features <- features
  }
  features <- features %||% rownames(x = object)
  
  if (is.null(cells.1) & !is.null(ident.1)){
    cells.1 <- Cells(subset(object, ident = ident.1))
  }
  if (is.null(cells.2) & !is.null(ident.2)){
    cells.2 <- Cells(subset(object, ident = ident.2))
  } else if (is.null(cells.2) & is.null(ident.2)){
    cells.2 <- Cells(subset(object, ident = ident.1, invert = T))
  }
  
  # error checking
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
  
  # feature selection (based on percentages)
  data <- switch(
    EXPR = slot,
    'scale.data' = counts,
    object
  )
  data <- data$ADT
  
  if (is.null(x = reduction)) {
    thresh.min <- 0
    pct.1 <- round(
      x = rowSums(x = data[features, cells.1] > thresh.min, na.rm = TRUE) /
        rowSums(!is.na(data[features, cells.1])),
      digits = 3
    )
    pct.2 <- round(
      x = rowSums(x = data[features, cells.2] > thresh.min, na.rm = TRUE) /
        rowSums(!is.na(data[features, cells.2])),
      digits = 3
    )
    data.alpha <- cbind(pct.1, pct.2)
    colnames(x = data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > min.pct))
    if (length(x = features) == 0) {
      stop("No features pass min.pct threshold")
    }
    alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
    features <- names(
      x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
    )
    if (length(x = features) == 0) {
      stop("No features pass min.diff.pct threshold")
    }
  } else {
    data.alpha <- data.frame(
      pct.1 = rep(x = NA, times = length(x = features)),
      pct.2 = rep(x = NA, times = length(x = features))
    )
  }
  
  # feature selection (based on average difference)
  mean.fxn <- if (is.null(x = reduction) && slot != "scale.data") {
    switch(
      EXPR = slot,
      'data' = function(x) {
        return(log(x = rowMeans(x = expm1(x = x), na.rm = TRUE) + pseudocount.use))
      },
      function(x) {
        return(log(x = rowMeans(x = x, na.rm = TRUE) + pseudocount.use))
      }
    )
  } else {
    rowMeans
  }
  data.1 <- mean.fxn(data[features, cells.1])
  data.2 <- mean.fxn(data[features, cells.2])
  
  total.diff <- (data.1 - data.2)
  
  if (is.null(x = reduction) && slot != "scale.data") {
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff > logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) > logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      stop("No features pass logfc.threshold threshold")
    }
  }
  
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    # Should be cells.1 and cells.2?
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }
  
  # perform DE
  if (!(test.use %in% c('negbinom', 'poisson', 'MAST', "LR")) && !is.null(x = latent.vars)) {
    warning(
      "'latent.vars' is only used for 'negbinom', 'poisson', 'LR', and 'MAST' tests",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  
  de.results <- switch(
    EXPR = test.use,
    'wilcox' = WilcoxDETest(
      data.use = object$ADT[features, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    stop("Unknown test: ", test.use)
  )
  
  if (is.null(x = reduction)) {
    diff.col <- ifelse(
      test = slot == "scale.data",
      yes = "avg_diff",
      no = "avg_logFC"
    )
    de.results[, diff.col] <- total.diff[rownames(x = de.results)]
    de.results <- cbind(de.results, data.alpha[rownames(x = de.results), , drop = FALSE])
  } else {
    diff.col <- "avg_diff"
    de.results[, diff.col] <- total.diff[rownames(x = de.results)]
  }
  if (only.pos) {
    de.results <- de.results[de.results[, diff.col] > 0, , drop = FALSE]
  }
  if (test.use == "roc") {
    de.results <- de.results[order(-de.results$power, -de.results[, diff.col]), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -de.results[, diff.col]), ]
    de.results$p_val_adj = p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = nrow(x = object)
    )
  }
  return(de.results)
}

# FindMarkers function for ADT
ADT_FindMarkers <- function(
  object,
  slot = "data",
  counts = numeric(),
  cells.1 = NULL, 
  cells.2 = NULL, 
  ident.1 = NULL,
  ident.2 = NULL,
  features = NULL,
  reduction = NULL,
  logfc.threshold = 0.15,
  test.use = 'wilcox',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  ...
) {
  de.results <- FindMarkers_na_removal(object = object,
                                       slot = slot,
                                       counts = counts,
                                       cells.1 = cells.1, 
                                       cells.2 = cells.2, 
                                       ident.1 = ident.1,
                                       ident.2 = ident.2,
                                       features = features,                                      
                                       reduction = reduction,
                                       logfc.threshold = logfc.threshold,
                                       test.use = test.use,
                                       min.pct = min.pct,
                                       min.diff.pct = min.diff.pct,
                                       verbose = verbose,
                                       only.pos = only.pos,
                                       max.cells.per.ident = max.cells.per.ident,
                                       random.seed = random.seed,
                                       latent.vars = latent.vars,
                                       min.cells.feature = min.cells.feature,
                                       min.cells.group = min.cells.group,
                                       pseudocount.use = pseudocount.use)
  
  de.results <- de.results[order(de.results$avg_logFC, decreasing=TRUE),]
  de.results <- de.results[de.results$p_val_adj < 0.05,]

  if (nrow(de.results) != 0){
    # AverageExpression
    tmp_avg_exp <- AverageExpression_na_removal(object, 
                                                features = rownames(de.results),
                                                assays = "ADT")$ADT
    colnames(tmp_avg_exp) <- paste(colnames(tmp_avg_exp), "_avg_exp", sep = "")
    de.results <- cbind(de.results, tmp_avg_exp)
    
    # row UMI counts
    if (!is.null(ident.1)){
      mtx1 <- object$ADT@counts[rownames(de.results), Cells(subset(object, ident = ident.1))] %>% as.matrix()
    } else{
      mtx1 <- object$ADT@counts[rownames(de.results), cells.1] %>% as.matrix()
    }
    if (!is.null(ident.2)){
      mtx2 <- object$ADT@counts[rownames(de.results), Cells(subset(object, ident = ident.2))] %>% as.matrix()
    } else{
      mtx2 <- object$ADT@counts[rownames(de.results), cells.2] %>% as.matrix()
    }
    
    if (ncol(mtx1) == 1){
      de.results[["ident.1_UMI_mean"]]   <- mean(mtx1, na.rm = T)
      de.results[["ident.1_UMI_median"]] <- median(mtx1, na.rm = T)
    } else{
      de.results[["ident.1_UMI_mean"]]   <- rowMeans(mtx1, na.rm = T)
      de.results[["ident.1_UMI_median"]] <- rowMedians(mtx1, na.rm = T)
    }
    if (ncol(mtx2) == 1){
      de.results[["ident.2_UMI_mean"]]   <- mean(mtx2, na.rm = T)
      de.results[["ident.2_UMI_median"]] <- median(mtx2, na.rm = T)
    } else{
      de.results[["ident.2_UMI_mean"]]   <- rowMeans(mtx2, na.rm = T)
      de.results[["ident.2_UMI_median"]] <- rowMedians(mtx2, na.rm = T)
    } 
    
    # add UMI counts and normalized counts data
    for (tmp_level in levels(object)){
      tmp_mtx <- object$ADT@data[rownames(de.results), Cells(subset(object, ident = tmp_level))] %>% as.matrix()  
      if (ncol(tmp_mtx) == 1){
        de.results[[paste(tmp_level, "_normalized_mean", sep = "")]] <- mean(tmp_mtx, na.rm = T)
        de.results[[paste(tmp_level, "_normalized_median", sep = "")]] <- median(tmp_mtx, na.rm = T)
        de.results[[paste(tmp_level, "_normalized_min", sep = "")]] <- min(tmp_mtx, na.rm = T)
        de.results[[paste(tmp_level, "_normalized_max", sep = "")]] <- max(tmp_mtx, na.rm = T)
      }else{
        de.results[[paste(tmp_level, "_normalized_mean", sep = "")]] <- rowMeans(tmp_mtx, na.rm = T)
        de.results[[paste(tmp_level, "_normalized_median", sep = "")]] <- rowMedians(tmp_mtx, na.rm = T)
        de.results[[paste(tmp_level, "_normalized_min", sep = "")]] <- rowMins(tmp_mtx, na.rm = T)
        de.results[[paste(tmp_level, "_normalized_max", sep = "")]] <- rowMaxs(tmp_mtx, na.rm = T)
      }
    }
    for (tmp_level in levels(object)){
      tmp_mtx2 <- object$ADT@counts[rownames(de.results), Cells(subset(object, ident = tmp_level))] %>% as.matrix()  
      if (ncol(tmp_mtx2) == 1){
        de.results[[paste(tmp_level, "_UMI_mean", sep = "")]] <- mean(tmp_mtx2, na.rm = T)
        de.results[[paste(tmp_level, "_UMI_median", sep = "")]] <- median(tmp_mtx2, na.rm = T)
        de.results[[paste(tmp_level, "_UMI_min", sep = "")]] <- min(tmp_mtx2, na.rm = T)
        de.results[[paste(tmp_level, "_UMI_max", sep = "")]] <- max(tmp_mtx2, na.rm = T)
      } else{
        de.results[[paste(tmp_level, "_UMI_mean", sep = "")]] <- rowMeans(tmp_mtx2, na.rm = T)
        de.results[[paste(tmp_level, "_UMI_median", sep = "")]] <- rowMedians(tmp_mtx2, na.rm = T)
        de.results[[paste(tmp_level, "_UMI_min", sep = "")]] <- rowMins(tmp_mtx2, na.rm = T)
        de.results[[paste(tmp_level, "_UMI_max", sep = "")]] <- rowMaxs(tmp_mtx2, na.rm = T)
      }
    }
  }
  # Merge
  de.results$avg_log2FC <- log2(exp(1)) * de.results$avg_logFC
  de.results$avg_log10FC <- log10(exp(1)) * de.results$avg_logFC
  de.results <- de.results[order(de.results$avg_logFC, decreasing=TRUE),]
  
  return(de.results)
}


# FindMarkers function for RNA
RNA_FindMarkers <- function (object,
                             ident.1 = NULL,
                             ident.2 = NULL,
                             group.by = NULL,
                             subset.ident = NULL,
                             assay = NULL,
                             slot = "data",
                             reduction = NULL,
                             features = NULL,
                             logfc.threshold = 0.25,
                             test.use = "wilcox",
                             min.pct = 0.1,
                             min.diff.pct = -Inf,
                             verbose = TRUE,
                             only.pos = TRUE,
                             max.cells.per.ident = Inf,
                             random.seed = 1,
                             latent.vars = NULL,
                             min.cells.feature = 3,
                             min.cells.group = 3,
                             pseudocount.use = 1){
  
  # FindMarkers
  tmp <- FindMarkers(object,
                     ident.1 = ident.1,
                     ident.2 = ident.2,
                     group.by = group.by,
                     subset.ident = subset.ident,
                     assay = assay,
                     slot = slot,
                     reduction = reduction,
                     features = features,
                     logfc.threshold = logfc.threshold,
                     test.use = test.use,
                     min.pct = min.pct,
                     min.diff.pct = min.diff.pct,
                     verbose = verbose,
                     only.pos = only.pos,
                     max.cells.per.ident = max.cells.per.ident,
                     random.seed = random.seed,
                     latent.vars = latent.vars,
                     min.cells.feature = min.cells.feature,
                     min.cells.group = min.cells.group,
                     pseudocount.use = pseudocount.use) 
  tmp$avg_log2FC <- log2(exp(1)) * tmp$avg_logFC
  tmp$avg_log10FC <- log10(exp(1)) * tmp$avg_logFC
  return(tmp)
}



RNA_FindAllMarkers <- function(
  object,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = 'wilcox',
  slot = 'data',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  return.thresh = 1e-2,
  ...
) {
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
    if (!PackageCheck('ape', error = FALSE)) {
      stop(cluster.ape, call. = FALSE)
    }
    tree <- Tool(object = object, slot = 'BuildClusterTree')
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
    descendants <- MapVals(
      vec = descendants,
      from = all.children,
      to = tree$tip.label
    )
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(
      node,
      as.numeric(x = setdiff(x = descendants, y = keep.children))
    )
    tree <- ape::drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }
  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        FindMarkers(
          object = object,
          assay = assay,
          ident.1 = if (is.null(x = node)) {
            idents.all[i]
          } else {
            tree
          },
          ident.2 = if (is.null(x = node)) {
            NULL
          } else {
            idents.all[i]
          },
          features = features,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          slot = slot,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          verbose = verbose,
          only.pos = only.pos,
          max.cells.per.ident = max.cells.per.ident,
          random.seed = random.seed,
          latent.vars = latent.vars,
          min.cells.feature = min.cells.feature,
          min.cells.group = min.cells.group,
          pseudocount.use = pseudocount.use,
          ...
        )
      },
      error = function(cond) {
        return(cond$message)
      }
    )
    if (is.character(x = genes.de[[i]])) {
      messages[[i]] <- genes.de[[i]]
      genes.de[[i]] <- NULL
    }
  }
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else if (is.null(x = node) || test.use %in% c('bimod', 't')) {
        gde <- gde[order(gde$p_val, -gde[, 2]), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = gde.all[, 2] > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  if (length(x = messages) > 0) {
    warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
    for (i in 1:length(x = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " versus all:\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }
  if (!is.null(x = node)) {
    gde.all$cluster <- MapVals(
      vec = gde.all$cluster,
      from = new.nodes,
      to = orig.nodes
    )
  }
  gde.all$avg_log2FC <- log2(exp(1)) * gde.all$avg_logFC
  gde.all$avg_log10FC <- log10(exp(1)) * gde.all$avg_logFC
  
  return(gde.all)
}




ADT_FindAllMarkers <- function(
  object,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = 'wilcox',
  slot = 'data',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  return.thresh = 1e-2,
  ...
) {
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
    if (!PackageCheck('ape', error = FALSE)) {
      stop(cluster.ape, call. = FALSE)
    }
    tree <- Tool(object = object, slot = 'BuildClusterTree')
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
    descendants <- MapVals(
      vec = descendants,
      from = all.children,
      to = tree$tip.label
    )
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(
      node,
      as.numeric(x = setdiff(x = descendants, y = keep.children))
    )
    tree <- ape::drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }
  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        FindMarkers_na_removal(
          object = object,
          assay = assay,
          cells.1 = NULL, 
          cells.2 = NULL,
          ident.1 = if (is.null(x = node)) {
            idents.all[i]
          } else {
            tree
          },
          ident.2 = if (is.null(x = node)) {
            NULL
          } else {
            idents.all[i]
          },
          features = features,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          slot = slot,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          verbose = verbose,
          only.pos = only.pos,
          max.cells.per.ident = max.cells.per.ident,
          random.seed = random.seed,
          latent.vars = latent.vars,
          min.cells.feature = min.cells.feature,
          min.cells.group = min.cells.group,
          pseudocount.use = pseudocount.use
        )
      },
      error = function(cond) {
        return(cond$message)
      }
    )
    if (is.character(x = genes.de[[i]])) {
      messages[[i]] <- genes.de[[i]]
      genes.de[[i]] <- NULL
    }
    if (nrow(genes.de[[i]]) != 0){
      # row UMI counts
      tmp_ident1 <- if (is.null(x = node)) {
        idents.all[i]
      } else {
        tree
      }
      tmp_ident2 <- if (is.null(x = node)) {
        NULL
      } else {
        idents.all[i]
      }
      if (!is.null(tmp_ident1)){
        mtx1 <- object$ADT@counts[rownames(genes.de[[i]]), Cells(subset(object, ident = tmp_ident1))] %>% as.matrix()
      } else{
        mtx1 <- object$ADT@counts[rownames(genes.de[[i]]), cells.1] %>% as.matrix()
      }
      if (!is.null(tmp_ident2)){
        mtx2 <- object$ADT@counts[rownames(genes.de[[i]]), Cells(subset(object, ident = tmp_ident2))] %>% as.matrix()
      } else{
        mtx2 <- object$ADT@counts[rownames(genes.de[[i]]), Cells(subset(object, ident = tmp_ident1, invert = T))] %>% as.matrix()
      }
      
      if (ncol(mtx1) == 1){
        genes.de[[i]][["ident.1_UMI_mean"]]   <- mean(mtx1, na.rm = T)
        genes.de[[i]][["ident.1_UMI_median"]] <- median(mtx1, na.rm = T)
      } else{
        genes.de[[i]][["ident.1_UMI_mean"]]   <- rowMeans(mtx1, na.rm = T)
        genes.de[[i]][["ident.1_UMI_median"]] <- rowMedians(mtx1, na.rm = T)
      }
      if (ncol(mtx2) == 1){
        genes.de[[i]][["ident.2_UMI_mean"]]   <- mean(mtx2, na.rm = T)
        genes.de[[i]][["ident.2_UMI_median"]] <- median(mtx2, na.rm = T)
      } else{
        genes.de[[i]][["ident.2_UMI_mean"]]   <- rowMeans(mtx2, na.rm = T)
        genes.de[[i]][["ident.2_UMI_median"]] <- rowMedians(mtx2, na.rm = T)
      } 
    }
  }
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else if (is.null(x = node) || test.use %in% c('bimod', 't')) {
        gde <- gde[order(gde$p_val, -gde[, 2]), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = gde.all[, 2] > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  if (length(x = messages) > 0) {
    warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
    for (i in 1:length(x = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " versus all:\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }
  if (!is.null(x = node)) {
    gde.all$cluster <- MapVals(
      vec = gde.all$cluster,
      from = new.nodes,
      to = orig.nodes
    )
  }
  gde.all$avg_log2FC <- log2(exp(1)) * gde.all$avg_logFC
  gde.all$avg_log10FC <- log10(exp(1)) * gde.all$avg_logFC
  
  return(gde.all)
}


CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots"),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message(
          "Suggested parameter: ",
          paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
          "\n"
        )
      }
    }
  }
}

FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  object <- UpdateSlots(object = object)
  slots <- na.omit(object = Filter(
    f = function(x) {
      sobj <- slot(object = object, name = x)
      return(is.list(x = sobj) && !is.data.frame(x = sobj) && !is.package_version(x = sobj))
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(inherits(x = object[[i]], what = classes.keep))
    }
  )
  object.classes <- which(x = object.classes, useNames = TRUE)
  return(names(x = object.classes))
}

UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(mode = class(x = xobj), length = 1L)
    }
  }
  return(object)
}

AverageExpression_na_removal <- function(
  object,
  assays = NULL,
  features = NULL,
  add.ident = NULL,
  slot = 'data',
  use.scale = FALSE,
  use.counts = FALSE,
  verbose = TRUE,
  ...
) {
  CheckDots(..., fxns = 'CreateSeuratObject')
  
  fxn.average <- switch(
    EXPR = slot,
    'data' = function(x) {
      rowMeans(x = expm1(x = x), na.rm = TRUE)
    },
    rowMeans
  )
  object.assays <- FilterObjects(object = object, classes.keep = 'Assay')
  assays <- assays %||% object.assays
  ident.orig <- Idents(object = object)
  orig.levels <- levels(x = Idents(object = object))
  ident.new <- c()
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(assays) == 0) {
      stop("None of the requested assays are present in the object")
    } else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (!is.null(x = add.ident)) {
    new.data <- FetchData(object = object, vars = add.ident)
    new.ident <- paste(
      Idents(object)[rownames(x = new.data)],
      new.data[, 1],
      sep = '_'
    )
    Idents(object, cells = rownames(new.data)) <- new.ident
  }
  data.return <- list()
  for (i in 1:length(x = assays)) {
    data.use <- GetAssayData(
      object = object,
      assay = assays[i],
      slot = slot
    )
    features.assay <- features
    if (length(x = intersect(x = features, y = rownames(x = data.use))) < 1 ) {
      features.assay <- rownames(x = data.use)
    }
    data.all <- list(data.frame(row.names = features.assay))
    for (j in levels(x = Idents(object))) {
      temp.cells <- WhichCells(object = object, idents = j)
      features.assay <- unique(x = intersect(x = features.assay, y = rownames(x = data.use)))
      if (length(x = temp.cells) == 1) {
        data.temp <- (data.use[features.assay, temp.cells])
        # transform data if needed (alternative: apply fxn.average to single value above)
        # if (!(use.scale | use.counts)) { # equivalent: slot.use == "data"
        if (slot == 'data') {
          data.temp <- expm1(x = data.temp)
        }
      }
      if (length(x = temp.cells) > 1 ) {
        data.temp <- fxn.average(data.use[features.assay, temp.cells, drop = FALSE])
      }
      data.all[[j]] <- data.temp
      if (verbose) {
        message(paste("Finished averaging", assays[i], "for cluster", j))
      }
      if (i == 1) {
        ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
      }
    }
    names(x = ident.new) <- levels(x = Idents(object))
    data.return[[i]] <- do.call(cbind, data.all)
    names(x = data.return)[i] <- assays[[i]]
  }
  return(data.return)
}
