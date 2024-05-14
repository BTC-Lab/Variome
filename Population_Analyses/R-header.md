

```r
# INIT --------------------------------------------------------------------

library(tidyverse)
library(igraph)
library(ggpubr)


# HELPER FUNCTIONS --------------------------------------------------------
```

Calculate variance and generate df with labels for plotting

@param val.obj the eigenval file that plink generates when computing PCA
@param vec.obj The eigenvec file that plink generates when computing PCA
@return A metrics data.frame containing PC IDs, variance, axis min/max


```r
getPlinkPCAMetrics <- function(val.obj, vec.obj) {
  
  t.var = sum(val.obj$V1)
  
  out.df <- data.frame(
    pc.id = paste("PC", seq_len(length(val.obj$V1)), sep = ""),
    var.explained = round(val.obj$V1 / t.var * 100, 2)) %>%
    dplyr::mutate(var.label = paste0(pc.id, ": ", var.explained, "%")) %>%
    dplyr::mutate(axis.min = (apply(vec.obj[,2:41], 2, function(col)
      min(col)))) %>%
    dplyr::mutate(axis.max = (apply(vec.obj[,2:41], 2, function(col)
      max(col))))
  
  rownames(out.df) <- out.df$pc.id
  
  return(out.df)
}
```

Calculate variance and generate df with labels for plotting from PRCOMP

@param prcomp.obj R base prcomp output object
@return A metrics data.frame containing PC IDs, variance, axis min/max


```r
getPCAMetrics <- function(prcomp.obj) {
  
  out.df <- data.frame(
    pc.id = paste("PC", seq_len(dim(prcomp.obj$x)[2]), sep = ""),
    var.explained =
      round((100 * prcomp.obj$sdev^2)/(sum(prcomp.obj$x^2/max(1,nrow(prcomp.obj$x) - 1))), 2)) %>%
    dplyr::mutate(var.label = paste0(pc.id, ": ", var.explained, "%")) %>%
    dplyr::mutate(axis.min = floor(apply(prcomp.obj$x, 2, function(col)
      min(col)))) %>%
    dplyr::mutate(axis.max = ceiling(apply(prcomp.obj$x, 2, function(col)
      max(col)))) %>%
    tibble::column_to_rownames(., var = "pc.id")
  
  return(out.df)
}
```

Add covariate information in meta.obj as attributes to the nodes of a graph

@param graph.obj The graph to add metadata to
@param meta.obj A data.frame of metadata, the rows must be named like the nodes
@return A graph with added information


```r
mergeMeta2Graph <- function(graph.obj, meta.obj) {
    
    for( c.idx in colnames(meta.obj) ) {
        graph.obj <- set_vertex_attr(graph.obj, name=eval(c.idx), value = meta.obj[V(graph.obj)$name, eval(c.idx)])
    }

    return(graph.obj)
}
```

Calculate metric for variability in a categorical variable

@param input.vec Vector of categorical variables, e.g., covariates or genetic variants
or a vector of relative frequencies
@param freq Boolean flag to indicate if input contains precomputed frequencies
@param standardize Boolean flag to indicate wether ouput should be standardized
@return A positive value indicating variability for given input. Higher is better.


```r
catVar <- function(input.vec, freq = FALSE, standardize = FALSE) {
    ## check if input is frequency table already
    if( !freq ) {
        input.vec <- table(input.vec) / length(input.vec)
        names(input.vec) <- paste("F", 1:length(input.vec), sep = "")
    }
    ## compute metric
    ssq_rf <- input.vec^2
    out.v <- 1 - sqrt(1 - (1- sum(ssq_rf)))
    ## check if standardization is required
    if( standardize ) {
        out.v <- out.v / (1 - (1 / sqrt(length(input.vec))))
    }
    return(out.v)
}
```

Find all clusters that contain a cycle and add an additional attribute to the
output of igraph::components that indicates if a cluster contains a cycle.
@params graph.obj the input graph to check for cycles in clusters
@return A list object containing cluster assignments and cycle indicator


```r
markCycle <- function( graph.obj, cutoff = -1 ) {
    cycle <- FALSE
    ncom <- igraph::components(graph.obj)
    ## init new cycle column with FALSE
    ncom$cycle <- rep(FALSE, times = ncom$no)
    ## for all clusters in the graph
    for( c.idx in 1:ncom$no ) {
        # get subgraph
        tmp.g <- igraph::induced.subgraph(graph.obj, vids = which(ncom$membership == c.idx))
        ## for all nodes in the subgraph
        for( v.idx in igraph::V(tmp.g) ) {
            ## find neighbors
            tmp.n <- igraph::neighbors(tmp.g, v.idx)
            ## for all neighbors / until break find paths
            for( n.idx in tmp.n ) {
                if( length(igraph::all_simple_paths(tmp.g, from = v.idx, to = n.idx, cutoff = eval(cutoff))) > 1 ) {
                    ## list longer than 1 means min. two patha have been found
                    ## two paths to the neighbor indicates a cycle
                    cycle <- TRUE
                    break
                }
            }
        ## if cycle we can stop searching
            if( cycle ) {
                ncom$cycle[c.idx] <- TRUE
                cycle <- FALSE
                break
            }
        }
        ## reach here and no cycle was found in this cluster
    }

    return(ncom)
}                                         
                                           
                                           
gStats <- function( graph.obj ) {
    tmp.ncom <- igraph::components(graph.obj)
    table(tmp.ncom$csize)
    #print("edge types")
    #table(E(graph.obj)$type)
}                                           
                                           

# PLOTTING FUNCTIONS ------------------------------------------------------

panelPCA <- function(plot.df, metric.df, p.axes, fill.var = NULL, col.var = NULL, l.pos = "none") {
  
  tmp.return <- ggplot(data=plot.df, aes(x=get(p.axes[1]), y=get(p.axes[2]), fill=get(fill.var), colour=get(col.var))) +
    geom_point(size=5, alpha = 0.75) + 
    scale_shape_manual(values=c(21, 22)) +
    theme_bw() +
    ## uncomment for admixture
    #scale_fill_manual(values = popCols, labels = c('African','American','East Asian','European','South Asian')) +
    #scale_colour_manual(values = popCols, name = "",
    #                    labels = c("Component 1","Component 2","Component 3","Component 4","Component 5"), 
    #                    breaks = c("European","South Asian","African","American","East Asian")) +
    xlim(metric.df[eval(p.axes[1]), "axis.min"], metric.df[eval(p.axes[1]), "axis.max"]) +
    ylim(metric.df[eval(p.axes[2]), "axis.min"], metric.df[eval(p.axes[2]), "axis.max"]) +
    xlab(metric.df[eval(p.axes[1]), "var.label"]) +
    ylab(metric.df[eval(p.axes[2]), "var.label"]) + 
    guides(shape="none",fill="none", colour = guide_legend(override.aes = list(size = 10))) +
    theme(legend.position=eval(l.pos), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x = element_text(size =26),
          axis.title.y = element_text(size =26),
          panel.border=element_blank(),
          legend.title =  element_text(size =26),
          legend.text =  element_text(size =26))
  
  return(tmp.return)
  
}
                                           
## quick and dirty graph plotting                                           
pG <- function(graph.obj) {
    
    plot.igraph( graph.obj,
              newpage=F, vertex.label.cex=0, vertex.label.color="#0057b8", vertex.label=V(graph.obj)$sampleID,
             edge.width=2.5,
              #edge.color = "#0057b8", #diffuStats::scores2colours(E(trio.graph)$weight, palette = colorRampPalette(c("blue", "white", "red"))),
             vertex.color="#b87500",
              mark.shape=1, mark.expand=20, vertex.size=1.5,
              main = paste0("some graph"),
              sub = "")

}
```

