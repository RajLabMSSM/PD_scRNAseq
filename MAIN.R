

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




ENSG_to_HGNC <- function(ensembl_ids, reference_genome="grch37", input="ensembl_gene_id"){ 
  print("+ Converting Ensembl IDs to HGNC Gene Symbols") 
  # Remove info after dot (numbers after dot specify the transcript)
  ensembl_ids <- gsub("\\..*","",ensembl_ids) 
  if(reference_genome=="grch38"){
    host <- "www.ensembl.org"
    } else {host <- paste0(reference_genome,".ensembl.org")}
  # BIOMART
  mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                          host=host, 
                          dataset = "hsapiens_gene_ensembl") 
  # View(biomaRt::listFilters(mart))
  # View(biomaRt::listAttributes(mart))
  res <- biomaRt::getBM(filters = input, # "ensembl_transcript_id" 
                        attributes = c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol"), 
                        values = ensembl_ids,
                        mart = mart)
  # Translate from feature to gene (or transcript)
  gene_dict <- setNames(res$hgnc_symbol, nm = res$ensembl_gene_id) 
  return(gene_dict)
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
                         results_path="./Results/DGE_results.txt",
                         force_DGE=F,
                         plot_topN=F,
                         plot_volcano=T){
  dge.start <- Sys.time()
  # Subset data
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
  
  # Import DGE results file if it already exists
  if(file.exists(results_path) & force_DGE==F){
    print("Results file already exists. Importing...")
    res <- data.table::fread(results_path)
  }else{ 
    # With regression:
    print("Initiating DGE analysis...")
    gene_fits <- monocle3::fit_models(cds_dge, 
                                      model_formula_str = paste0("~",variable), 
                                      expression_family = "quasipoisson",
                                      cores = nCores, 
                                      verbose = T)
    fit_coefs <- coefficient_table(gene_fits) 
    res <- fit_coefs %>% filter(term != "(Intercept)") %>% arrange(q_value, desc(abs(estimate)))  
    # Save
    if(results_path!=F){
      data.table::fwrite(res, results_path)
    }
    # Report
    ## Time taken
    dge.end <- Sys.time()
    dge.diff <- as.numeric(round(difftime(dge.end, dge.start, units = "hours"),2))
    print(paste("+ Calculated in",dge.diff,"hours."))
    ## Sig genes
    sig_genes <- res %>% filter (q_value < 0.05) %>% pull(gene_short_name)
    print(paste("DGE:",dim(cds_DGE)[1],"genes x",dim(cds_DGE)[2],"samples"))
    print(paste("+",length(sig_genes),"significiant DGE(s) detected.")) 
  } 
  
  # Plot
  if(plot_topN!=F){
    try({
      top_genes <- res$gene_short_name[1:plot_topN] %>%  as.character()
      p <- monocle3::plot_genes_violin(cds_dge[top_genes,], 
                                       group_cells_by=variable, 
                                       ncol=2)
      print(p)
    }) 
  } 
  
  if(plot_volcano){
    try({
      varis <- unique(pData(cds_dge)[variable])[,1]
      comparison_label <- paste0(variable," : ", varis[1], " vs. ", varis[2] )
      vp <- volcano_plot(dge = res, caption=comparison_label, topN_labeled=6)
      print(vp)
    }) 
  } 
  return(res)
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
  # Y max
  neg.log.pvals <- -log10(dge$q_value)
  yMaximum <- max(neg.log.pvals)
  yMax  <- yMaximum + yMaximum/8 #ifelse(max(-log10(dge$p_val_adj))<45, 50, max(-log10(dge$p_val_adj)) + 10)
  # X max
  dge_sig <- subset(dge, q_value<0.05)
  xMaximum <- max(abs(dge_sig$estimate))
  xMax <- xMaximum + xMaximum/8
  
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
    xlim(-xMax, xMax) +
    ylim(0, yMax) + 
    theme(plot.title = element_text(hjust = 0.5))
  return(vol)
}



top_cluster_markers <- function(cds_DGE, 
                                cluster_list=c(1,2), 
                                genes_to_test_per_group=100, 
                                save_path= "./Results/cluster_markers.csv", 
                                verbose=F){
  # Subset clusters
  cds_clusts <- cds_DGE[,(pData(cds_DGE)$Cluster %in% cluster_list)]
  # Test
  marker_test_res = monocle3::top_markers(cds_clusts,
                                          group_cells_by="cluster", 
                                          genes_to_test_per_group = genes_to_test_per_group,
                                          cores=nCores, 
                                          verbose=verbose)
  # Sort
  marker_test_res <- marker_test_res %>% arrange(desc(fraction_expressing), marker_test_q_value, desc(specificity), desc(pseudo_R2))  
  # Save
  if(save_path!=F){
    data.table::fwrite(marker_test_res, save_path)
  }
  # Filter and plot top markers
  ## Filter
  filtered_res <- marker_test_res %>%
    filter(fraction_expressing >= 0.10 & marker_test_q_value <=0.05) %>%
    group_by(cell_group) %>% arrange(desc(specificity), desc(pseudo_R2)) 
  ## Plot
  top_specific_markers <- filtered_res$gene_id[1:min(c(10,nrow(filtered_res)))] %>% unique()
  print(paste("Plotting the top",length(top_specific_markers), "specific markers."))
  monocle3::plot_cells(cds_DGE, 
                       genes = top_specific_markers, 
                       show_trajectory_graph = F)
  # Return full results
  return(marker_test_res)
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



zscore <- function(vec, rescale_01=T){ 
  vec.sd <- sd(vec)*sqrt((length(vec)-1)/(length(vec)))
  vec.mean <- mean(vec)
  z <- (vec - vec.mean) / vec.sd
  if(rescale_01){
    z <- scales::rescale(x=z, to=c(1,0))
  }
  return(z)
}

run_enrichr <- function(DEG_table, 
                        results_path = "./Results/topEnrichR_hits.csv",  
                        enrichr_dbs = c("KEGG_2018", 
                                        "Reactome_2016",
                                        "GO_Biological_Process_2018", 
                                        "GO_Molecular_Function_2018",
                                        "GO_Cellular_Component_2018",
                                        "Rare_Diseases_AutoRIF_ARCHS4_Predictions",
                                        "ARCHS4_Cell-lines",
                                        "Aging_Perturbations_from_GEO_up", 
                                        "Aging_Perturbations_from_GEO_down",
                                        "Human_Gene_Atlas",
                                        "ChEA_2016", 
                                        "KEA_2015"),
                        top_genes=500,
                        include_weights=T
){ 
  # createDT(enrichR::listEnrichrDbs(), "Enrichr Databases") 
  DEGs_table <- subset(DEGs_table, q_value<=0.05)[1:top_genes,]
  if(include_weights){
    vec <- -log10(DEGs_table$q_value)
    vec[is.infinite(vec)] <- .Machine$double.xmin
    z <- zscore(vec, rescale_01 = T) 
    geneList <- data.frame(Gene=row.names(DEGs_table), 
                           Weight=z)
  } else {
    geneList <- data.frame(Gene=row.names(DEGs_table), 
                           Weight=1)
  }
 
  results <- enrichR::enrichr(genes = geneList, databases = enrichr_dbs) 
  
  topHits <- data.table()
  for (db in enrichr_dbs){
    cat('\n')
    cat("##",db,"\n")  
    # res <- subset(results[[db]], Adjusted.P.value<=0.05) 
    newdf = subset(results[[ db ]], Adjusted.P.value<=0.05, select=c("Term","Overlap","Adjusted.P.value")) %>% mutate(Database=db)  
    newdf = newdf[1:3, ]
    topHits <- rbind(topHits, newdf[complete.cases(newdf), ]) 
    cat('\n')
  }  
  write.csv(topHits, results_path)
  return(topHits)
}



report_overlap <- function(genomeSize, list1, list2, verbose=T){  
  go.obj <-  GeneOverlap::newGeneOverlap(listA = list1, listB = list2, genome.size = genomeSize)
  go.obj <-  GeneOverlap::testGeneOverlap(go.obj) 
  
  overlappingGenes <-  GeneOverlap::getIntersection(go.obj)
  percent_of_targetGenes <- length(overlappingGenes) / length(list2)*100
  percent_of_DEGs <- length(overlappingGenes) / length(list1)*100
  targetGenes_DEGs <- list2[overlappingGenes %in% list2]
  if(verbose){
    print(go.obj)
    cat("\n",round(percent_of_targetGenes, 2),"% of list2 genes (",length(overlappingGenes),"/",length(list2),") are in list1.")
    
    cat("\n",round(percent_of_DEGs, 2),"% of the list1 genes (",
        length(overlappingGenes),"/",length(list1),") are in list2.")
    
    cat("\n-------------------------------------------------------",
        "\n\n++++++++++++ Enrichment p-value =",go.obj@pval,"++++++++++++\n\n")
  } 
  return(list(go.obj=go.obj,
              overlappingGenes=overlappingGenes )) 
} 


module_vs_module_enrichment <- function(mod.df1, 
                                        mod.df2, 
                                        genomeSize, 
                                        verbose=T,
                                        save_file="./Results/bulk.modules_vs_sc.modules.txt"
                                        ){ 
  mod1.unique <- length(unique(mod.df1$module))
  mod2.unique <- length(unique(mod.df2$module))
  print(paste("mod.df1 contains", mod1.unique, "unique modules."))
  print(paste("mod.df2 contains", mod2.unique, "unique modules."))
  print(paste("Conducting enrichment tests on",mod1.unique*mod2.unique,"module-module combinations..."))
  start.mod <- Sys.time() 
  # Level 1
  MOD.dt <-lapply(unique(mod.df1$module), function(MOD){
    sub.df1 <- subset(mod.df1, module==MOD)
    # Level 2
    mod.dt <- lapply(unique(mod.df2$module), function(mod, MOD.=MOD){
      if(verbose){print(paste(MOD.,"(mod.df1) vs.",mod,"(mod.df2)"))}
      sub.df2 <- subset(mod.df2, module==mod)
       
      res.list <- report_overlap(genomeSize = genomeSize,
                                 list1 = sub.df1$gene.symbol, 
                                 list2 = sub.df2$gene.symbol,  
                                 verbose = F)
      go.obj <- res.list$go.obj
      # overlappingGenes <- res.list$overlappingGenes
      summary.dt <- data.table::data.table(Module1.name = MOD, 
                                           Module2.name = mod, 
                                           Module1.size = length(go.obj@listA), 
                                           Module2.size = length(go.obj@listB), 
                                           Module1.proportion.overlap = round(length(go.obj@intersection) / length(go.obj@listA),3),
                                           Module2.proportion.overlap = round(length(go.obj@intersection) / length(go.obj@listB),3),
                                           intersection.size = length(go.obj@intersection), 
                                           union.size = length(go.obj@union),
                                           genome.size = go.obj@genome.size,
                                           p.value = go.obj@pval, 
                                           odds.ratio = go.obj@odds.ratio,
                                           Jaccard.index = go.obj@Jaccard)
    }) %>% data.table::rbindlist()
  }) %>% data.table::rbindlist() 
  
  MOD.dt <- MOD.dt %>% dplyr::mutate(FDR=p.adjust(p.value, method="fdr"), 
                                     Bonferroni=p.adjust(p.value, method="bonferroni")) %>% arrange(FDR)
  MOD.sig <- MOD.dt %>% subset(FDR<=0.05) 
  end.mod <- Sys.time()
  print(paste(nrow(MOD.dt),"enrichment tests conducted in", round(end.mod-start.mod,2),"seconds."))
  print(paste(nrow(MOD.sig),"enrichment tests were significant (at FDR ≤ 0.05)."))
  # MOdule conservation
  mod.df1.convervation <- length(unique(MOD.sig$Module1.name)) / length(unique(MOD.dt$Module1.name))
  mod.df2.convervation <- length(unique(MOD.sig$Module2.name)) / length(unique(MOD.dt$Module2.name))
  print(paste0(round(mod.df1.convervation*100,2),"% of mod.df1 modules showed enrichment for a mod.df2 module."))
  print(paste0(round(mod.df2.convervation*100,2),"% of mod.df2 modules showed enrichment for a mod.df1 module."))
  
  if(save_file!=F){
    data.table::fwrite(MOD.dt, save_file)
  }
  
  return(MOD.dt)
}


module_vs_module_enrichment.terms <- function(MOD.dt, 
                                              mod.enrich, 
                                              top.terms=F,
                                              row.limit=F){
  print(paste("Comparing overlap of enrichment terms for all module-module combinations."))
  # Set limits 
  if(row.limit==F){row.limit <- nrow(MOD.dt) }
  
  
  mm.overlap <- lapply(1:row.limit, function(i){
    MOD.dt <- (MOD.dt %>% arrange(desc(Jaccard.index))) 
    
    # Module 1
    mod1 <- MOD.dt[i,"Module1.name"]
    mod1.terms <- (subset(mod.enrich, module==mod1) %>%  
                     arrange(FDR, desc(precision), desc(recall)))
    if(top.terms==F){top.terms.m1 <- nrow(mod1.terms) }else{top.terms.m1 <- top.terms}
    mod1.terms <- mod1.terms[1:top.terms.m1,]$term_name %>% unique()
    # Module 2
    mod2 <- MOD.dt[i,"Module2.name"]
    mod2.terms <- (subset(mod.enrich, module==mod2) %>%  
                     arrange(FDR, desc(precision), desc(recall)))
    if(top.terms==F){top.terms.m2 <- nrow(mod2.terms) }else{top.terms.m2 <- top.terms}
    mod2.terms <- mod2.terms[1:top.terms.m2,]$term_name %>% unique()
    
    # Test overlap
    union.terms <- union_all(mod1.terms, mod2.terms)
    go.obj <-  GeneOverlap::newGeneOverlap(listA = mod1.terms, 
                                           listB = mod2.terms, 
                                           genome.size = length(union.terms))
    go.obj <-  GeneOverlap::testGeneOverlap(go.obj) 
    
    summary.dt <- data.table::data.table(Module1.name = mod1, 
                                         Module2.name = mod2, 
                                         Module1.size = length(go.obj@listA),
                                         Module2.size = length(go.obj@listB), 
                                         Module1.top.term = go.obj@listA[1],
                                         Module2.top.term = go.obj@listB[1],
                                         Module1.proportion.overlap = round(length(go.obj@intersection) / length(go.obj@listA),3),
                                         Module2.proportion.overlap = round(length(go.obj@intersection) / length(go.obj@listB),3),
                                         intersection.size = length(go.obj@intersection), 
                                         union.size = length(go.obj@union),
                                         ontology.size = go.obj@genome.size,
                                         p.value = go.obj@pval, 
                                         odds.ratio = go.obj@odds.ratio,
                                         Jaccard.index = go.obj@Jaccard)
    return(summary.dt)
  }) %>% data.table::rbindlist(fill=T)
  return(mm.overlap)
}



gprofiler2.module_enrichment <- function(mod.df, dataset=""){
  start.enrich <- Sys.time()
  mod.enrich <- lapply(unique(mod.df$module), function(mod){ 
    print(paste("gprofiler2:: Running enrichment on module:",mod))
    gostres <- gprofiler2::gost(query = subset(mod.df, module==mod)$gene.symbol, organism = "hsapiens") 
    if(is.null(gostres)){
      enrich.res <- data.frame(dataset=dataset, module=mod)
    } else{
      enrich.res <- gostres$result
      enrich.res <- cbind(dataset=dataset, module=mod, enrich.res)
    } 
    return(enrich.res)
  }) %>% data.table::rbindlist(fill = T)
  
  end.enrich <- Sys.time()
  print(paste("gprofiler2::",length(unique(mod.df$module)),
              "tested for ontological enrichment in",round(end.enrich-start.enrich,2),"seconds"))
  return(mod.enrich)
}


# HEATMAP
overlap_heatmap <- function(clustDAT, geneList, 
                            title="Overlapping Genes:\nGene List vs. DGE Genes", 
                            interactive=F){ 
  knitr::opts_chunk$set(fig.width=14, fig.height=14) 
  markerDF  <- get_markerDT(clustDAT, markerList = geneList, rawData = F)
  markerMatrix <- reshape2::acast(markerDF, Gene~Cluster, value.var="Expression",
                                  fun.aggregate = mean, drop = F, fill = 0)
  # Heatmap.2 
  my_palette <- colorRampPalette(c("purple", "black", "cyan"))(n = 1000)
  hmap <- gplots::heatmap.2(markerMatrix, xlab = "Cluster", dendrogram = "row",
                            col = my_palette, tracecol = "gray", srtCol = 0, offsetCol=1.5, vline=T, 
                            trace='none', key.title=NA, key.ylab = "Expression", colsep=T, sepwidth = 0.01, 
                            main = title) 
  knitr::opts_chunk$set(fig.width=7, fig.height=5) 
}

overlap_expression_plot <- function(DEG_table, geneList, title=""){
  dat.sub <- subset(DEG_table, gene_short_name %in% geneList) %>% arrange(desc(estimate))
  dat.sub$gene_short_name <- factor(dat.sub$gene_short_name, levels = dat.sub$gene_short_name)
  
  p <- ggplot(dat.sub, aes(x=gene_short_name, y=estimate, fill=estimate)) + 
    geom_col(show.legend = F) + labs(title=title) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) #+  ylim(-2,2)
  plot(p)
}




# ^^^^^^^^^^^^ 3D Volcano Plot ^^^^^^^^^^^^ 
volcano_3d <- function(dge_path="./Results/across_Clust1.vs.Clust2.csv"){
  library(plotly)
  dge = data.table::fread("./Results/across_Clust1.vs.Clust2.csv")
  norm <- function(x) {(x - min(x)) / (max(x) - min(x))}
  # dge$estimate.norm <-  norm(dge$estimate)
  # dge$q_value.norm <-  norm(-log10(dge$q_value))
  # dge$colorscale <- dge$estimate.norm + dge$q_value.norm
  # dge[is.na(dge)] <-0
  # dge$cone <-  pi*dge$estimate#((dge$estimate)^2)*(-log10(dge$q_value)/3) 
  # dge$cone.surface = pi*dge$estimate*(dge$estimate+sqrt(dge$q_value^2+dge$estimate^2))
  p <- plot_ly(dge, 
               x = ~test_val, 
               y = ~estimate, 
               z = ~-log10(q_value),
               text = ~ gene_short_name,
               marker = list( color = colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(200), showscale = T)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'Normalized Effect'),
                        yaxis = list(title = 'Effect'),
                        zaxis = list(title = '-log(q-value)'))) 
  p 
  
}





