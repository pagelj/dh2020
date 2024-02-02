#!/usr/bin/env Rscript

library(DramaAnalysis)
library(tidyverse)
library(igraph)
library(ggplot2)

source("configuration.R") # tz6r has NA values in segmentColumn, using a modified version of the configuration function of DramaAnalysis

ids <- grep("gdc-coref", loadAllInstalledIds(), value = TRUE)
ids <- ids[! ids %in% c("gdc-coref:nks0.0")] # Remove Cato, because text is empty
dramas <- loadDrama(ids)
dramas <- split(dramas)

plot_degree_curve <- function(dramas, segment = "Scene", n = 3) {
  lapply(dramas,
         function(d) {
           c_active <- configuration(d, onlyPresence = TRUE, mode = "Active", segment = segment) %>%
             filterCharacters(d, n = n) %>%
             characterNames(d)
           c_passive <- configuration(d, onlyPresence = TRUE, mode = "Passive", segment = segment) %>%
             filterCharacters(d, n = n) %>%
             characterNames(d)
           segs_active <- 4:ncol(c_active)
           segs_passive <- 4:ncol(c_passive)
           degree_active <- lapply(segs_active,
                                   function(s) {
                                     cs <- c_active[c(1:3,s)]
                                     mat <- as.matrix(cs)
                                     copr <- mat %*% t(mat)
                                     rownames(copr) <- c_active$character
                                     colnames(copr) <- c_active$character
                                     g <- igraph::graph_from_adjacency_matrix(copr, 
                                                                              weighted=TRUE,
                                                                              mode="undirected",
                                                                              diag=FALSE
                                     )
                                     degree <- igraph::strength(g)
                                     degree
                                   }
           )
           degree_passive <- lapply(segs_passive,
                                    function(s) {
                                      cs <- c_passive[c(1:3,s)]
                                      mat <- as.matrix(cs)
                                      copr <- mat %*% t(mat)
                                      rownames(copr) <- c_passive$character
                                      colnames(copr) <- c_passive$character
                                      g <- igraph::graph_from_adjacency_matrix(copr, 
                                                                               weighted=TRUE,
                                                                               mode="directed",
                                                                               diag=FALSE
                                      )
                                      degree <- igraph::strength(g)
                                      degree
                                    }
           )
           df_active <- as.data.frame(do.call(cbind,degree_active))
           df_active$mode <- "Co-Presence"
           df_active <- tibble::rownames_to_column(df_active, var = "character")
           df_passive <- as.data.frame(do.call(cbind,degree_passive))
           df_passive$mode <- "Co-Reference"
           df_passive <- tibble::rownames_to_column(df_passive, var = "character")
           df <- rbind(df_active, df_passive)
           df <- gather(df, segment, value, 2:(ncol(df)-1))
           df$segment <- gsub("V", "", df$segment)
           df$segment <- as.integer(df$segment)
           plot_name <- dramaNames(d)
           gg <- ggplot(data = df) + 
             aes(x = segment, y = value, group = character) +
             facet_grid(mode ~ character) +
             geom_line() +
             ylim(c(0, max(df$value))) +
             ggtitle(plot_name) +
             xlab(segment) +
             ylab("Weighted Degree") +
             theme_bw()
           ggsave(plot = gg, 
                  filename = paste0("plots/", "centrality_curve", "_", plot_name, "_", "degree", ".png"),
                  dpi = 300)
           gg
         }
  )
}
plots <- plot_degree_curve(dramas = dramas, segment = "Scene", n = 3)

get_centrality_measures <- function(dramas, segment = "Scene", mode = "Active") {
  centrality <- lapply(dramas,
                       function(d) {
                         c <- configuration(d, onlyPresence = TRUE, mode = mode, segment = segment) %>%
                           characterNames(d)
                         mat <- as.matrix(c)
                         copr <- mat %*% t(mat)
                         rownames(copr) <- c$character
                         colnames(copr) <- c$character
                         g <- igraph::graph_from_adjacency_matrix(copr, 
                                                                  weighted=TRUE,
                                                                  mode="undirected",
                                                                  diag=FALSE
                         )
                         degree <- igraph::degree(g, normalized = TRUE)
                         #strength <- igraph::strength(g)
                         betweenness <- igraph::betweenness(g, normalized = TRUE)
                         closeness <- igraph::closeness(g, normalized = TRUE)
                         eigen_centrality <- igraph::eigen_centrality(g)$vector
                         centrality <- rbind(degree, betweenness, closeness, eigen_centrality)
                         centrality
                       }
  )
  df <- as.data.frame(t(as.data.frame(centrality)))
  df
}

centralities_active <- get_centrality_measures(dramas, segment = "Scene", mode = "Active")
centralities_passive <- get_centrality_measures(dramas, segment = "Scene", mode = "Passive")
active <- summarise(centralities_active, 
          mean.degree = mean(degree), 
          mean.betweenness = mean(betweenness), 
          mean.closeness = mean(closeness), 
          mean.eigen_centrality = mean(eigen_centrality),
          sd.degree = sd(degree, na.rm = TRUE), 
          sd.betweenness = sd(betweenness, na.rm = TRUE), 
          sd.closeness = sd(closeness, na.rm = TRUE), 
          sd.eigen_centrality = sd(eigen_centrality, na.rm = TRUE))
active$mode <- "co-presence"
passive <- summarise(centralities_passive, 
          mean.degree = mean(degree), 
          mean.betweenness = mean(betweenness), 
          mean.closeness = mean(closeness), 
          mean.eigen_centrality = mean(eigen_centrality),
          sd.degree = sd(degree), 
          sd.betweenness = sd(betweenness), 
          sd.closeness = sd(closeness), 
          sd.eigen_centrality = sd(eigen_centrality))
passive$mode <- "co-reference"
table <- rbind(active, passive)
table <- pivot_longer(table, cols = 1:8, names_to = "centrality", values_to = "value")
table[c('method', 'centrality')] <- str_split_fixed(table$centrality, '\\.', 2)
table <- pivot_wider(table, names_from = "method")
table