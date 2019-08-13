

# Render Datatables
createDT <- function(DF, caption="", scrollY=500){
  data <- DT::datatable(DF, caption=caption,
                        extensions =  'Buttons',
                        options = list( dom = 'Bfrtip', 
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,  
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  ) 
  return(data)
}
## Need to wrap inside tagList when rendering DTs within a for loop
createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY)) 
}

# Import parameters supplied through yaml header
import_parameters <- function(params){    
  params_list <- list(
    resultsPath=params$resultsPath,
    # nCores=as.numeric(params$nCores),
    subsetGenes=params$subsetGenes,
    subsetCells=params$subsetCells,
    resolution=params$resolution,
    perplexity=as.numeric(params$perplexity)
  )
  # Have to setwd via knitr 
  knitr::opts_chunk$set(echo=T, error=T, root.dir = params$resultsPath 
                        # cache=T, cache.lazy=T
  )  
  # Utilize parallel processing later on
  cat("**** __Utilized Cores__ **** =", params$nCores)
  
  bindingIsLocked("params", env = .GlobalEnv)
  unlockBinding("params", env = .GlobalEnv)
  print(params)  
  for (n in names(params_list) ){
    assign(n, params_list[n],  env = .GlobalEnv)
  }
}


save_3D_clusters <- function(cds, output_path="./3D_objects"){
  dir.create(output_path, showWarnings = F, recursive = T)
  for(clust in levels(clusters(cds)) ){ 
    print(paste0("Creating OBJ for cluster ",clust))
    cds_sub <- cds[,clusters(cds)==clust] 
    rgl::plot3d(x = cds_sub@reducedDims$UMAP[,1],
                y = cds_sub@reducedDims$UMAP[,2],
                z = cds_sub@reducedDims$UMAP[,3],  
                box=F, axes=F)
    rgl::writeOBJ(con = paste0(output_path,"/cluster",clust,".obj"),  )
  } 
}


seurat_to_monocle <- function(seurat_object, monocle_version=3){ 
  # Make sure Seurat object is up to do date with Seurat v3 format
  sDAT <- Seurat::UpdateSeuratObject(object = DAT)
   
  
  # From: https://github.com/cole-trapnell-lab/monocle-release/issues/262
  print("Processing...")
  print("+ Expression data")
  expr_dat <- sDAT@raw.data
  
  print("+ Phenotype data") 
  # Filter metadata 
  sDAT@meta.data  <- sDAT@meta.data[!is.na(sDAT@meta.data$nGene),]  
  metadata <- sDAT@meta.data
  
  print("+ Feature Data")
  fDATA <- data.frame(gene_short_name = row.names(expr_dat), row.names = row.names(expr_dat))
  
  
  # Make sure metadata and expression data have the same cells
  cell_IDs <- dplyr::intersect(colnames(expr_dat), row.names(metadata))
  expr_dat <- expr_dat[,cell_IDs]
  metadata <- metadata[cell_IDs,]
  print("Expression Data dims:")
  print(dim(expr_dat))
  print("Metadata dims:")
  print(dim(metadata))
  print("Feature Data dims:")
  print(dim(fDATA))
  
  
  #Construct monocle cds  
  if(monocle_version==3){
    print("+ Converting to monocle (Version 3)")
    cds <- monocle3::new_cell_data_set(expression_data = expr_dat,
                                               cell_metadata = metadata,
                                               gene_metadata = fDATA)
  } else {
    print("+ Converting to monocle (Version 2)")
    pd <- new('AnnotatedDataFrame', data = metadata) 
    fd <- new('AnnotatedDataFrame', data = fDATA)
    cds <- monocle::newCellDataSet(cellData = expr_dat,
                                   phenoData = pd,
                                   featureData = fd,
                                   lowerDetectionLimit = 0.5,
                                   expressionFamily = VGAM::negbinomial.size())
}
  
  return(cds)
}



get_biotypes <- function(gene_list, subsetGenes="protein_coding"){
  # gene_list <- row.names(DAT@raw.data) 
    cat(paste("Subsetting genes:",subsetGenes, "\n"))
    # If the gene_biotypes file exists, import csv. Otherwise, get from biomaRt
    if(file_test("-f", file.path(root, "Data/gene_biotypes.csv"))){
      biotypes <- read.csv(file.path(root, "Data/gene_biotypes.csv"))
    }
    else {
      ensembl <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                                  dataset="hsapiens_gene_ensembl") 
      # ensembl <- biomaRt::useDataset(mart = ensembl, dataset = "hsapiens_gene_ensembl")
      biomaRt::listFilters(ensembl)
      biomaRt::listAttributes(ensembl)   
      biotypes <- biomaRt::getBM(attributes=c("hgnc_symbol", "gene_biotype"), 
                                 filters="hgnc_symbol",
                                 values=gene_list, 
                                 mart=ensembl) 
      write.csv(biotypes, file.path(root,"Data/gene_biotypes.csv"), quote=F, row.names=F)
    } 
  return(biotypes)
}


subset_seurat <- function(DAT, genes.use){
  subset.matrix <- DAT@raw.data[genes.use, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
  subDAT <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
  orig.ident <- DAT@ident # Pull the identities from the original Seurat object as a data.frame
  subDAT <- AddMetaData(object = subDAT, metadata = orig.ident, col.name = "meta.data") # Add the idents to the meta.data slot
  subDAT <- SetAllIdent(object = subDAT, id = "ident") # Assign identities for the new Seurat object
  return(subDAT)
}

monocle3_DGE <- function(cds_DGE, 
                         variable="dx", 
                         nCores=4, 
                         expression_family="quasipoisson", 
                         variable_subsets=F,
                         plot_topN=F,
                         plot_volcano=T,
                         results_path="./Results/DGE_results.txt"){
  # With regression:
  suppressWarnings(
  if(variable_subsets==F){
    cds_dge <- cds_DGE
  } else { 
    # variable="Cluster"
    # variable_subsets = c(1,2) 
    bool_vector <- (pData(cds_DGE)[variable] %in% variable_subsets)[[1]]
    cds_dge <- cds_DGE[, bool_vector]
  }
  )
  dge.start <- Sys.time()
  gene_fits <- fit_models(cds_dge, 
                          model_formula_str = paste0("~",variable), 
                          expression_family = "quasipoisson",
                          cores = nCores, 
                          verbose = T)
  fit_coefs <- coefficient_table(gene_fits) 
  res <- fit_coefs %>% filter(term != "(Intercept)") %>% arrange(q_value, desc(abs(estimate)))
  # time_terms <- time_terms %>% mutate(q_value = p.adjust(p_value))
  sig_genes <- res %>% filter (q_value < 0.05) %>% pull(gene_short_name)
  
  # Report
  dge.end <- Sys.time()
  dge.diff <- as.numeric(round(difftime(dge.end, dge.start, units = "hours"),2))
  print(paste("DGE:",dim(cds_DGE)[1],"genes x",dim(cds_DGE)[2],"samples"))
  print(paste("+ Calculated in",dge.diff,"hours."))
  print(paste("+",length(sig_genes),"significiant DGE(s) detected."))
  
  # Plot
  if(plot_topN!=F){
    top_genes <- res$gene_short_name[1:plot_topN] %>%  as.character()
    p <- monocle3::plot_genes_violin(cds_dge[top_genes,], 
                                     group_cells_by=variable, 
                                     ncol=2)
    print(p)
  } 
    
  if(plot_volcano){
    varis <- unique(pData(cds_dge)[variable])[,1]
    comparison_label <- paste0(variable," : ", varis[1], " vs. ", varis[2] )
    vp <- volcano_plot(dge = res, caption=comparison_label, topN_labeled=6)
    print(vp)
  }
  
  if(save_results!=F){
    data.table::fwrite(paste0("./Results/",results_path))
  }
 
  return(res)
  # # With graph autocorrelation:
  # pr_test_res <- graph_test(cds[variance[1:10],],  
  #                           neighbor_graph="principal_graph", 
  #                           cores=nCores)
  # pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
}


volcano_plot <- function(dge, caption="", topN_labeled=6){
  # dge <- res
  # Replace zeros with smallest representable number in R to prevent turning into Infs
  dge[dge$q_value==0, "q_value"] <- .Machine$double.xmin
  dge$sig<-  ifelse( dge$q_value<0.05 & abs(dge$estimate)<=1.5, "q_value<0.05",
                     ifelse( dge$q_value<0.05  & abs(dge$estimate)>=1.5, "q_value<0.05 & estimate>=1.5",
                             "q_value>0.05"
                     ))  
  dge <- arrange(dge, desc(sig))
  neg.log.pvals <- -log10(dge$q_value)
  maximum <- max(neg.log.pvals)
  yMax  <- maximum + maximum/8 #ifelse(max(-log10(dge$p_val_adj))<45, 50, max(-log10(dge$p_val_adj)) + 10)
  
  vol <- ggplot(data=dge, aes(x=estimate, y= -log10(q_value))) +
    geom_point(alpha=0.5, size=3, aes(col=sig)) + 
    scale_color_manual(values=list("q_value<0.05"="turquoise3",
                                   "q_value<0.05 & estimate>=1.5"="purple", 
                                   "q_value>0.05" = "darkgray")) +
    theme(legend.position = "none") + 
    theme_bw() + 
    xlab(expression(paste("Estimate"))) +
    ylab(expression(paste(-log^{10},"(q_value)"))) +  
    ## ggrepl labels
    geom_text_repel(data=arrange(dge,  q_value, desc(estimate))[1:topN_labeled,], 
                    # filter(dge, avg_logFC>=1.5)[1:10,],
                    aes(label=gene_short_name),  color="black", alpha=.5,
                    segment.color="black", segment.alpha=.5, seed = 2019) +  
    # Lines
    geom_vline(xintercept= -1.5,lty=4, lwd=.3, alpha=.5) + 
    geom_vline(xintercept= 1.5,lty=4, lwd=.3, alpha=.5) +
    geom_hline(yintercept= -log10(0.05),lty=4, lwd=.3, alpha=.5) +   
    ggtitle(caption) +
    xlim(-2,2) + ylim(0, yMax) + 
    theme(plot.title = element_text(hjust = 0.5))
  return(vol)
}

get_earliest_principal_node <- function(cds,  variable = "dx", variable_value = "PD"){
  cell_ids <- which(colData(cds)[, variable] == variable_value)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


