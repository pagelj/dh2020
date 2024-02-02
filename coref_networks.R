#!/usr/bin/env Rscript

library(DramaAnalysis)
library(tidyverse)
library(igraph)
library(ggplot2)

source("configuration.R") # tz6r has NA values in segmentColumn, using a modified version of the configuration function of DramaAnalysis

output_dir <- "plots/"

ids <- grep("gdc-coref", loadAllInstalledIds(), value = TRUE)
#ids <- ids[! ids %in% c("gdc-coref:nks0.0")] # Remove Cato, because text is empty
dramas <- loadDrama(ids)
dramas <- split(dramas)

plot_networks <- function(dramas, segment = "Scene", mode = "Active", save = FALSE) {
  graphs <- lapply(seq(dramas),
                   function(i) {
                     d <- dramas[[i]]
                     c <- configuration(d, onlyPresence = TRUE, mode = mode, segment = segment) %>%
                       characterNames(d)
                     mat <- as.matrix(c)
                     copr <- mat %*% t(mat)
                     rownames(copr) <- c$character
                     colnames(copr) <- c$character
                     g <- igraph::graph_from_adjacency_matrix(copr, 
                                                              weighted=TRUE,
                                                              mode=ifelse(mode == "Active", "undirected", "directed"),
                                                              diag=FALSE
                     )
                     plot_name <- dramaNames(d)
                     if (save) {
                       png(filename = paste0(output_dir, "network_", plot_name, "_", mode, "_", segment, ".png"),
                           width = 8, height = 8, units = 'in', res = 300)
                     }
                     plot(g,
                          layout = layout.circle(g),
                          edge.label = E(g)$weight,
                          edge.width = 1,
                          edge.arrow.size = 0.5,
                          vertex.size = scales::rescale(to = c(2, 15), igraph::strength(g, mode = "all")),
                          vertex.label.dist=1,
                          vertex.color = adjustcolor("SkyBlue2", alpha.f = .5),
                          vertex.label.color = adjustcolor("black", alpha.f = .8),
                          edge.curved=ifelse(mode == "Active", 0.0, 0.4),
                          main = plot_name
                     )
                     if (save) {
                       dev.off()
                     }
                     g
                   }
  )
  graphs
}
plots_active <- plot_networks(dramas, segment = "Scene", mode = "Active", save = TRUE)
plots_passive <- plot_networks(dramas, segment = "Scene", mode = "Passive", save = TRUE)