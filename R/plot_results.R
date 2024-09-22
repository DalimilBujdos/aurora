plot_results <- function(pheno_mat, dist_mat, color, df_fianl_rf) {
  color <- data.frame(pheno = levels(pheno_mat$pheno),
                      colors = c("red", "blue", "green", "black", "yellow"))
  df_final_rf <- read_csv("C:/Users/bujdo/OneDrive/Dokumenty/PhD_OToole/main_project/data/reuteri/allo_auto_multiclass/aurora_results_rf_2022-11-27.csv")
  df_final_rf <- as.data.frame(df_final_rf)
  # size is proportional the the classfication median
  # each phenotype should have its own color
  min_max_norm <- function(x, min, max) {
    # function performing min/max normalization
    x <- 1 + ( ((x - min)*(5-1))/ (max - min) )
    return(x)
  }

  dist_mat <- final_proxy_mat_sum
  mds <- cmdscale(dist_mat, eig = TRUE)
  mds_df <- as.data.frame(mds$points)
  mds_df$keys <- rownames(mds_df)
  # find the right color for each strain
  mds_df$pheno <- pheno_mat$pheno
  mds_df$color <- color$colors[match(mds_df$pheno, color$pheno)]
  # find the right size for each strain
  mds_df$size <- rep(NA, nrow(mds_df))
  for (i in 1:nrow(mds_df)) {
    strain <- mds_df$keys[i]
    strain_pheno <- mds_df$pheno[i]
    strain_pheno <- lookup_tbl[as.numeric(strain_pheno)]
    search_lbl <- paste0(strain_pheno, "_classification_median")
    mds_df$size[i] <- df_final_rf[df_final_rf$strain_ID == strain,
                                  colnames(df_final_rf) == search_lbl]
  }
  mds_df$size <- as.numeric(mds_df$size)
  mds_df$size <- abs(mds_df$size - 1)
  mds_df$size <- sapply(mds_df$size, min_max_norm, min(mds_df$size), max(mds_df$size))

  colnames(mds_df)[c(1,2)] <- c("dim_1", "dim_2")

  point_plot <- ggplot() +
    geom_point(data = mds_df, aes(x = dim_1, y = dim_2), color = mds_df$color, size = mds_df$size) +
    scale_colour_manual(values=setNames(color$colors, color$pheno)) +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    theme(axis.text = element_text(size = 20)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(legend.position="none")
  plot(point_plot)

  # plot heatmap

  heatmap(dist_mat, ColSideColors = mds_df$color, RowSideColors = mds_df$color, margins = c(1,1), keep.dendro = TRUE)

  ##############################################################################
  # plot results from adaboost
  min_max_norm <- function(x, min, max) {
    # function performing min/max normalization
    x <- 1 + ( ((x - min)*(5-1))/ (max - min) )
    return(x)
  }
  color <- data.frame(pheno = levels(pheno_mat$pheno),
                      colors = c("red", "blue", "green", "black", "yellow"))
  df_final_ada <- read_csv("C:/Users/bujdo/OneDrive/Dokumenty/PhD_OToole/main_project/data/reuteri/allo_auto_multiclass/aurora_results_ada_2022-11-27.csv")
  df_final_ada <- as.data.frame(df_final_ada)
  # get all the classification medians into one df
  index <- seq(from = 2, to = ncol(df_final_ada), by = 3)
  df_dist <-  df_final_ada[, index]
  # calculate manhattan dist
  dist_mat <- as.matrix(dist(as.matrix(df_dist[,-ncol(df_dist)]),
                        method = "minkowski",
                        diag = TRUE, upper = TRUE))

  mds <- cmdscale(dist_mat, eig = TRUE)
  mds_df <- as.data.frame(mds$points)
  mds_df$keys <- df_final_ada$strain_ID
  # find the right color for each strain
  mds_df$pheno <- pheno_mat$pheno
  mds_df$color <- color$colors[match(mds_df$pheno, color$pheno)]
  # find the right size for each strain
  mds_df$size <- rep(NA, nrow(mds_df))
  for (i in 1:nrow(mds_df)) {
    strain <- mds_df$keys[i]
    strain_pheno <- mds_df$pheno[i]
    strain_pheno <- lookup_tbl[strain_pheno]
    search_lbl <- paste0(strain_pheno, "_classification_median")
    mds_df$size[i] <- df_final_ada[df_final_ada$strain_ID == strain,
                                  colnames(df_final_ada) == search_lbl]
  }
  mds_df$size <- as.numeric(mds_df$size)
  mds_df$size <- abs(mds_df$size - 1)
  mds_df$size <- sapply(mds_df$size, min_max_norm, min(mds_df$size), max(mds_df$size))

  colnames(mds_df)[c(1,2)] <- c("dim_1", "dim_2")

  point_plot <- ggplot() +
    geom_point(data = mds_df, aes(x = dim_1, y = dim_2), color = mds_df$color, size = mds_df$size) +
    scale_colour_manual(values=setNames(color$colors, color$pheno)) +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    theme(axis.text = element_text(size = 20)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(legend.position="none")
  plot(point_plot)

  heatmap(dist_mat, ColSideColors = mds_df$color, RowSideColors = mds_df$color, margins = c(1,1), keep.dendro = TRUE)

  ##############################################################################
  # plot misslabel results from aura and rf
  # plot results from adaboost
  min_max_norm <- function(x, min, max) {
    # function performing min/max normalization
    x <- 1 + ( ((x - min)*(5-1))/ (max - min) )
    return(x)
  }
  colnames(df_threshold_rf) <- c(as.character(phenotypes), "pheno")
  df_threshold_rf$pheno <- as.factor(df_threshold_rf$pheno)
  color <- data.frame(pheno = 1:length(phenotypes),
                      colors = c("red", "blue", "green", "pink", "yellow"))

  dist_mat <- as.matrix(dist(as.matrix(df_threshold_rf[,-ncol(df_threshold_rf)]),
                             method = "manhattan",
                             diag = TRUE, upper = TRUE))

  mds <- cmdscale(dist_mat, eig = TRUE)
  mds_df <- as.data.frame(mds$points)
  # find the right color for each strain
  mds_df$pheno <- df_threshold_rf$pheno
  mds_df$color <- rep(NA, nrow(mds_df))
  mds_df$color <- color$colors[match(mds_df$pheno, color$pheno)]
  mds_df$color[is.na(mds_df$color)] <- "black"
  colnames(mds_df)[c(1,2)] <- c("dim_1", "dim_2")

  point_plot <- ggplot() +
    geom_point(data = mds_df, aes(x = dim_1, y = dim_2), color = mds_df$color) +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    theme(axis.text = element_text(size = 20)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(legend.position="none")
  plot(point_plot)

  heatmap(dist_mat, ColSideColors = mds_df$color, RowSideColors = mds_df$color, margins = c(1,1), keep.dendro = TRUE)


  return(list(point_plot = point_plot, heat_map = heat_map))
}
