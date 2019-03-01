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

# Raise Memory Limit 
## Rstudio has a default memory limit of only 1GB. To override this, detect the true memory available and set a new limit.
raise_memory_limit <- function(){ 
  # library(ulimit) # devtools::install_github("krlmlr/ulimit")
  # library(benchmarkme) 
  RAM <- print(benchmarkme::get_ram())
  ## Convert GB to Mib
  RAM_Mib <- strsplit(RAM, " ")[[1]][1] %>% as.numeric() * 953.67431640625
  cat(paste("Available RAM:",RAM))
  ## Set new memory limit 
  ulimit::memory_limit(RAM_Mib)  
}


subsetBiotypes <- function(DAT, subsetGenes="protein_coding"){
  if( subsetGenes!=F ){
    cat(paste("Subsetting genes:",subsetGenes, "\n"))
    # If the gene_biotypes file exists, import csv. Otherwise, get from biomaRt
    if(file_test("-f", file.path(root, "Data/gene_biotypes.csv"))){
      biotypes <- read.csv(file.path(root, "Data/gene_biotypes.csv"))
    }
    else {
      ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                         dataset="hsapiens_gene_ensembl") 
      ensembl <- useDataset(mart = ensembl, dataset = "hsapiens_gene_ensembl")
      listFilters(ensembl)
      listAttributes(ensembl)   
      biotypes <- getBM(attributes=c("hgnc_symbol", "gene_biotype"), filters="hgnc_symbol",
                        values=row.names(DAT@data), mart=ensembl) 
      write.csv(biotypes, file.path(root,"Data/gene_biotypes.csv"), quote=F, row.names=F)
    } 
    # Subset data by creating new Seurat object (annoying but necessary)
    geneSubset <- biotypes[biotypes$gene_biotype==subsetGenes,"hgnc_symbol"] 
    
    cat(paste(dim(DAT@raw.data[geneSubset, ])[1],"/", dim(DAT@raw.data)[1], 
              "genes are", subsetGenes))
    # Add back into DAT 
    subset.matrix <- DAT@raw.data[geneSubset, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
    DAT_sub <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
    orig.ident <- row.names(DAT@meta.data) # Pull the identities from the original Seurat object as a data.frame
    DAT_sub <- AddMetaData(object = DAT_sub, metadata = DAT@meta.data) # Add the idents to the meta.data slot
    DAT_sub <- SetAllIdent(object = DAT_sub, id = "ident") # Assign identities for the new Seurat object
    DAT <- DAT_sub
    rm(list = c("DAT_sub","geneSubset", "subset.matrix", "orig.ident")) 
  } 
  return(DAT)
}

remove_nonmatched_metadata <- function(DAT, subsetCells){
  # Get rid of any NAs (cells that don't match up with the metadata) 
  if(subsetCells==F){
    DAT <- FilterCells(object = DAT,  subset.names = "nGene", low.thresholds = 0)
  } else {DAT <- FilterCells(object = DAT,  subset.names = "nGene", low.thresholds = 0,
                             # Subset for testing 
                             cells.use = DAT@cell.names[0:as.numeric(subsetCells) ]
    )
  }  
  return(DAT)
}



 
volcanoPlot <- function(DEG_df, caption="", topFC_labeled=5){
  DEG_df$sig<-  ifelse( DEG_df$p_val_adj<0.05 & DEG_df$avg_logFC<1.5, "p_val_adj<0.05",
                        ifelse( DEG_df$p_val_adj<0.05  & DEG_df$avg_logFC>1.5, "p_val_adj<0.05 & avg_logFC>1.5",
                                "p_val_adj>0.05"
                        )) 
  DEG_df <- arrange(DEG_df, desc(sig))
  neg.log.pvals <- -log10(DEG_df$p_val_adj) 
  maximum <- max(neg.log.pvals[!is.infinite(neg.log.pvals)])
  yMax  <- maximum + maximum/8 #ifelse(max(-log10(DEG_df$p_val_adj))<45, 50, max(-log10(DEG_df$p_val_adj)) + 10)
  
  vol <- ggplot(data=DEG_df, aes(x=avg_logFC, y= -log10(p_val_adj))) +
    geom_point(alpha=0.5, size=3, aes(col=sig)) + 
    scale_color_manual(values=list("p_val_adj<0.05"="turquoise3",
                                   "p_val_adj<0.05 & avg_logFC>1.5"="purple", 
                                   "p_val_adj>0.05" = "darkgray")) +
    theme(legend.position = "none") + 
    xlab(expression(paste("Average ",log^{2},"(fold change)"))) +
    ylab(expression(paste(-log^{10},"(p-value)"))) + 
    xlim(-2,2) + ylim(0, yMax) +
    ## ggrepl labels
    geom_text_repel(data= arrange(DEG_df,  p_val_adj, desc(avg_logFC))[1:topFC_labeled,], 
                    # filter(DEG_df, avg_logFC>=1.5)[1:10,],
                    aes(label=gene),  color="black", alpha=.5,
                    segment.color="black", segment.alpha=.5  
    ) +  
    # Lines
    geom_vline(xintercept= -1.5,lty=4, lwd=.3, alpha=.5) + 
    geom_vline(xintercept= 1.5,lty=4, lwd=.3, alpha=.5) +
    geom_hline(yintercept= -log10(0.05),lty=4, lwd=.3, alpha=.5) + 
    ggtitle(caption) 
  print(vol)
}
 

get_markerDT <- function(DAT, markerList, rawData=T, meta_vars =c("barcode", "dx", "mut","Cluster", "percent.mito","nGene", "nUMI")){
  if(rawData==T){
    exp <- DAT@raw.data 
  }else{exp <- DAT@scale.data} 
  marker.matrix <- exp[row.names(exp) %in% markerList, ] 
  markerMelt <- reshape2:::melt.matrix(marker.matrix, varnames = c("Gene", "Cell"), value.name = "Expression")  
  markerMeta <- DAT@meta.data[meta_vars] %>% copy()  
  markerDT <- base::merge(markerMelt,markerMeta, by.x="Cell", by.y="barcode")  
  return(markerDT)
}

replace_zero_pvals <- function(DEGs){
  # convert 0 to the smallest number R can represent(2.225074e-308)
  DEGs[DEGs$p_val==0.0000000000000000000000000000,] <- 2.225074e-308  
  DEGs[DEGs$p_val_adj==0.0000000000000000000000000000,] <- 2.225074e-308
  return(DEGs)
}




runDGE <- function(DAT, meta_var, group1, group2, test.use="wilcox", show.table=T, show.volcano=T, allGenes=F){
  #print(paste("DGE_allCells",meta_var,sep="_")) 
  DAT <- SetAllIdent(DAT, id = meta_var)
  DAT <- StashIdent(DAT, save.name = meta_var)  
  if(allGenes==T){
    # DON'T use filtering and apply FDR
    DEGs <- FindMarkers(DAT, ident.1=group1, ident.2=group2, test.use=test.use, only.pos = F, 
                        logfc.threshold = 0, min.pct = 0, min.cells.group = 1, 
                        min.cells.gene = 1, min.diff.pct = -Inf)
  } else {
    # Use filtering and apply Bonferonni-correct p-vals
    DEGs <- FindMarkers(DAT, ident.1=group1, ident.2=group2, test.use=test.use, only.pos = F)
  } 
  DEGs <- replace_zero_pvals(DEGs)
  DEGs$gene <- row.names(DEGs)
  
  cap <- paste("DEGs:\n",group1, "vs.", group2)
  if(show.table==T){
    createDT_html(DEGs, caption = cap)
  }else{print("Not showing table")}
  if(show.volcano==T){ volcanoPlot(DEGs, caption = cap)
  }else(print("Not showing volcano plot")) 
  
  DAT <- SetAllIdent(DAT, id = "Cluster")
  return(DEGs)
}

DGE_within_clusters <- function(DAT, meta_var, group1, group2, clusterList, allClusts=F, allGenes = F){ 
  if(allClusts==T){ clusterList <- unique(DAT@meta.data$Cluster) }
  for (clust in clusterList){ 
    cat('\n')   
    cat("### ",paste("Cluster ",clust,": ",group1," vs. ", group2, sep="") , "\n")
    DAT_clustSub <- Seurat::SubsetData(DAT, subset.name ="Cluster", accept.value = clust, subset.raw = T)  
    DEG_df <-runDGE(DAT_clustSub, meta_var, group1, group2, allGenes) 
    cat('\n')   
  } 
}
