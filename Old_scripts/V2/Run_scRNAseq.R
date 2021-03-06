# Also try 'argparse' package for more complex handling

args = commandArgs(trailingOnly = T)
 
subsetGenes <- ifelse(args[1]==F, 'all', args[1])
subsetCells = args[2]
resolution = args[3] 
perplexity = args[4]
# detectCores() Identifies all cores on a cluster, even if you don't have access to them all. 
# This causes an error that forces multithread-enabled functions to revert to the default
nCores = args[5] #parallel::detectCores() 

#nCores <-args[4] #ifelse(is.null(args[4]), parallel::detectCores(), args[4]) 

parameters <- paste(sep='', 'subsetGenes-',subsetGenes,'__subsetCells-',subsetCells, 
                    '__Resolution-', resolution,'__perplexity-',perplexity,
                    '__nCores-',nCores)

resultsPath <- file.path('Results',parameters)
dir.create(resultsPath, showWarnings=F, recursive=T) 
dir.create(file.path("Data"), showWarnings=F) 
params_list <- list(subsetGenes= subsetGenes, subsetCells=subsetCells, resolution=resolution, 
                    resultsPath=resultsPath, nCores=nCores, perplexity=perplexity)
 
# 1) Preprocessing
# rmarkdown::render(input = 'scRNAseq_Preprocessing.Rmd',
#                   params = params_list,
#                   output_file = file.path(resultsPath, paste("Preprocessing.html",sep=""))
#                   )
# #


# 3) DGE
rmarkdown::render(input = 'scRNAseq_DGE.Rmd',
                  params = params_list,
                  output_file = file.path(resultsPath, paste("DGE.html",sep=""))
)

# # 4) About page
# rmarkdown::render(input = 'index.Rmd', 
#                   params = params_list,  
#                   output_file = file.path(resultsPath, paste("index.html",sep=""))
# )
