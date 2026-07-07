# Compatibility shims for Seurat v5 vs ggplot2 >=3.4 and Seurat S7 issues
# DOC: Seurat v5 uses S7 classes for BPCells/arrow assays; older plotting code uses S3 dispatch

if (packageVersion("ggplot2") >= "3.4.0") {
  local({
    original_facet_grid <- ggplot2::facet_grid
    new_facet_grid <- function(facets = NULL, rows = NULL, ...) {
      if (!is.null(facets) && is.null(rows)) rows <- facets
      original_facet_grid(rows = rows, ...)
    }
    assignInNamespace("facet_grid", new_facet_grid, ns = "ggplot2")
  })
  cat("Applied facet_grid compatibility patch for ggplot2 >=3.4\n")
}

# Safe wrappers - all visualization functions wrapped in tryCatch
safe_heatmap <- function(...) {
  tryCatch(DoHeatmap(...), error = function(e) {
    cat("DoHeatmap skipped:", conditionMessage(e), "\n")
    ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("DoHeatmap skipped")
  })
}

safe_dotplot <- function(...) {
  tryCatch(DotPlot(...), error = function(e) {
    cat("DotPlot skipped:", conditionMessage(e), "\n")
    ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("DotPlot skipped")
  })
}

safe_featureplot <- function(...) {
  tryCatch(FeaturePlot(...), error = function(e) {
    cat("FeaturePlot skipped:", conditionMessage(e), "\n")
    ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("FeaturePlot skipped")
  })
}
