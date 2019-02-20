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
  print( htmltools::tagList( createDT(DF, caption, scrollY)) ) 
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


volcanoPlot <- function(DEG_df, caption="", topFC_labeled=5){
  DEG_df$sig<-  ifelse( DEG_df$p_val_adj<0.05 & DEG_df$avg_logFC<1.5, "p_val_adj<0.05",
                        ifelse( DEG_df$p_val_adj<0.05  & DEG_df$avg_logFC>1.5, "p_val_adj<0.05 & avg_logFC>1.5",
                                "p_val_adj>0.05"
                        )) 
  DEG_df <- arrange(DEG_df, desc(sig))
  yMax  <- max(-log10(DEG_df$p_val_adj)) + max(-log10(DEG_df$p_val_adj))/3 #ifelse(max(-log10(DEG_df$p_val_adj))<45, 50, max(-log10(DEG_df$p_val_adj)) + 10)
  
  vol <- ggplot(data=DEG_df, aes(x=avg_logFC, y= -log10(p_val_adj))) +
    geom_point(alpha=0.5, size=3, aes(col=sig)) + 
    scale_color_manual(values=list("p_val_adj<0.05"="turquoise3",
                                   "p_val_adj<0.05 & avg_logFC>1.5"="purple", 
                                   "p_val_adj>0.05" = "darkgray")) +
    theme(legend.position = "none") + 
    xlab(expression(paste("Average ",log^{2},"(fold change)"))) +
    ylab(expression(paste(-log^{10},"(p-value)"))) + xlim(-2,2) + ylim(0, yMax) +
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
 

get_markerDT <- function(DAT, markerList, rawData=T, meta_vars =c("barcode", "dx", "mut","post_clustering", "percent.mito","nGene", "nUMI")){
  if(rawData==T){
    exp <- DAT@raw.data 
  }else{exp <- DAT@scale.data} 
  marker.matrix <- exp[row.names(exp) %in% markerList, ] 
  markerMelt <- reshape2:::melt.matrix(marker.matrix, varnames = c("Gene", "Cell"), value.name = "Expression")  
  markerMeta <- DAT@meta.data[meta_vars] %>% copy()  
  markerDT <- base::merge(markerMelt,markerMeta, by.x="Cell", by.y="barcode")  
  return(markerDT)
}
