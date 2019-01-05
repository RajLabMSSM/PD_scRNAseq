args = commandArgs(trailingOnly = T)

subsetGenes = args[1]
subGenes <- ifelse(subsetGenes==F, 'all', subsetGenes)
subsetCells = args[2]
resolution = args[3] 
#nCores <-args[4] #ifelse(is.null(args[4]), parallel::detectCores(), args[4]) 

parameters <- paste(sep='', 'subsetGenes-',subGenes,'__subsetCells-',subsetCells, '__Resolution-', resolution)
root = getwd()
resultsPath <- file.path(root,'Results',parameters)
dir.create(resultsPath, showWarnings=F, recursive=T)
        

rmarkdown::render(input = 'run_seurat.Rmd', 
                  params = list(subsetGenes= subsetGenes, subsetCells=subsetCells, resolution=resolution, resultsPath=resultsPath), 
                  output_file = file.path(resultsPath,paste(parameters,".html",sep=""))
                  )