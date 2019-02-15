# Import parameters supplied through yaml header
import_parameters <- function(params){    
  params_list <- list(
  resultsPath=params$resultsPath,
     nCores=params$nCores,
     subsetGenes=params$subsetGenes,
     subsetCells=params$subsetCells,
     resolution=params$resolution,
     perplexity=params$perplexity
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
