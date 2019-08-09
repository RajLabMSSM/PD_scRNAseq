
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
    pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data) 
    fd <- new('AnnotatedDataFrame', data = fDATA)
    cds <- monocle::newCellDataSet(data,
                                           phenoData = pd,
                                           featureData = fd,
                                           lowerDetectionLimit = 0.5,
                                           expressionFamily = VGAM::negbinomial.size())
  }
  
  return(cds)
}
