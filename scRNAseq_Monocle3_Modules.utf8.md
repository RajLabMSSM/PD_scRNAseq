---
title: "Gene Co-expression Modules"
subtitle: "Bulk vs. Single-Cell RNA-seq"
author: "<h3>Author</h3>Brian M. Schilder, Bioinformatician II"
date: "2019-10-09"
output: 
  rmarkdown::html_document: 
    theme: spacelab
    highlight: zenburn 
    code_folding: show 
    toc: true 
    toc_float: true
    smooth_scroll: true
    number_sections: false 
    self_contained: true ### !IMPORTANT!: Must == true in order to create GitHub page.
params:
  subsetGenes: "protein_coding" #FALSE
  subsetCells: 500 #FALSE
  resolution: 0.6
  resultsPath: "./Results"
  nCores: 2
  perplexity: 30
editor_options: 
  chunk_output_type: inline
---


```r
library(dplyr) 
library(ggplot2) 
root = "./"
source(file.path(root,"MAIN.R")) 
import_parameters(params)
```

```
## **** __Utilized Cores__ **** = 2$subsetGenes
## [1] "protein_coding"
## 
## $subsetCells
## [1] 500
## 
## $resolution
## [1] 0.6
## 
## $resultsPath
## [1] "./Results"
## 
## $nCores
## [1] 2
## 
## $perplexity
## [1] 30
```

```r
load("./Data/monocle3_CDS.RData") 
nCores <- 4
```

# Graph-autocorrelation Analysis 

Find genes of interest through graph-autocorrelation analysis.


```r
# pr_graph_test_res = monocle3::graph_test(cds, neighbor_graph="knn", cores=nCores)
# pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))
```

# Louvain Modules

- Identify co-expression modules.<br>
`find_gene_modules` essentially runs UMAP on the genes (as opposed to the cells) and then groups them into modules using Louvain community analysis.
- **Advantages**: Much faster than other alogrithms like WGCNA.  
- **Disadvantages**: Seems to be limited of ~1000 genes.


```r
gene_module_df = monocle3::find_gene_modules(cds[var.genes[1:1000], ], 
                                             # resolution=1e-2, 
                                             reduction_method = "UMAP",
                                             max_components = 3,
                                             cores=nCores,
                                             umap.fast_sgd = F, # Slower but determinstic results
                                             random_seed = 2019)
data.table::fwrite(gene_module_df, "./Results/sc.modules_by_gene.txt", sep="\t")
createDT(gene_module_df)
```

preserved4d917e2b9b06f54

## Pheatmap

- Plot the eigengene expression of each gene co-expression module (y-axis) within each cell cluster (x-axis).


```r
# Create pheatmap table
cell_group_df = tibble::tibble(cell=row.names(SummarizedExperiment::colData(cds)),
                               cell_group=monocle3::clusters(cds)[colnames(cds)])
agg_mat = monocle3::aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) = stringr::str_c("Cluster ", colnames(agg_mat))
module_cluster_exp <- data.frame(as.matrix(agg_mat))
data.table::fwrite(data.table::data.table(module_cluster_exp, keep.rownames = "Module"),
                   "./Results/sc.modules_by_cluster.txt", sep="\t" )
# Plot
pheatmap::pheatmap(agg_mat, cluster_rows=T, cluster_cols=F,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6) 
```

<img src="scRNAseq_Monocle3_Modules_files/figure-html/Gene Modules - pheatmap-1.png" width="672" />
 


## 3D Scatter


```r
# plotly R documentation: https://plot.ly/r/
# rgl::plot3d(x = gene_module_df$dim_1, y=gene_module_df$dim_2, z=gene_module_df$dim_3, ) 
p.3d <- plotly::plot_ly(gene_module_df, 
                     x = ~dim_1, y = ~dim_2, z = ~dim_3, 
                     color = ~module,  
                     text = ~paste("Gene:",id),
                     colors = "Dark2",  
                     # alpha = .7, 
                     size = 5
                     ) %>%
  plotly::layout(title = 'scRNA-seq Gene Co-expression Modules')  %>%
  plotly::add_markers() %>%
  plotly::layout(scene = list(xaxis = list(title = 'UMAP 1'),
                     yaxis = list(title = 'UMAP 2'),
                     zaxis = list(title = 'UMAP 3')))
htmltools::tagList(setNames(list(p.3d), NULL)) 
```

preserve0206a201d2f8518f


## Module Expression


```r
# gene_module_df <-data.table::fread("./Results/sc.modules_by_gene.txt")
monocle3::plot_cells(cds,
                     genes=gene_module_df,
                     group_cells_by="cluster",
                     color_cells_by="cluster",
                     show_trajectory_graph=F) #  rasterize = T (for super high-res figures)
```

<img src="scRNAseq_Monocle3_Modules_files/figure-html/Module Expression-1.png" width="672" />

```r
p2 <- monocle3::plot_cells(cds,
                           # sc.module 9 correponds to sc.module 5 in the old analyses
                           genes=subset(gene_module_df, module==9),
                           group_cells_by="cluster",
                           color_cells_by="cluster",
                           show_trajectory_graph=F)
print(p2)
```

<img src="scRNAseq_Monocle3_Modules_files/figure-html/Module Expression-2.png" width="672" />

```r
ggsave("./Results/sc.module9.png", plot=p2, dpi=400, width=5, height=3.5)
```





# Bulk RNAseq Modules vs. scRNAseq Modules

- Co-expression modules were identified through the [WGCNA analysis of bulk monocyte RNA-seq data](https://github.com/RajLabMSSM/WGCNA_MyND/tree/master/2nd_batch_Analysis/protein_coding/04_sva/deepSplit_4/signed), conducte by Katia Lopes.  
- Here we compare the modules identified in the scRNA-seq data here (**sc.modules**), to the modules identified in the bulk monocyte data (**bulk.modules**).  


```r
# Prepare bulk modules
bulk.modules <- data.table::fread("./Data/Katia_Modules/Gene_modules.tbl.txt") 
# bulk.modules$ensembl.id <- gsub("\\.*","", bulk.modules$gene_ids)
# ensembl.dict <- ENSG_to_HGNC(ensembl_ids = unique(bulk.modules$gene_ids), 
#                              reference_genome="grch38", 
#                              input = "ensembl_transcript_id")
if(!file.exists("./Data/Katia_Modules/ens.geneid.gencode.v30")){
  system("scp schilb03@chimera.hpc.mssm.edu:/sc/hydra/projects/pd-omics/GRCH38.p12/ens.geneid.gencode.v30 ./Data/Katia_Modules")
} 
ensembl.dt <- data.table::fread("./Data/Katia_Modules/ens.geneid.gencode.v30")
ensembl.dict <- setNames(ensembl.dt$GeneSymbol, ensembl.dt$gene_id)
bulk.modules$gene.symbol <- ensembl.dict[bulk.modules$gene_ids]  
bulk.modules <- bulk.modules %>% dplyr::rename(module=moduleColors)
createDT(bulk.modules)
```

```
## Warning in instance$preRenderHook(instance): It seems your data is too
## big for client-side DataTables. You may consider server-side processing:
## https://rstudio.github.io/DT/server.html
```

preserved5a04596609fbf5b

```r
# Prepare sc modules
sc.modules <- data.table::fread("./Results/sc.modules_by_gene.txt") %>% dplyr::rename(gene.symbol=id)
createDT(sc.modules)
```

preservedcd1a05bc218ba3e

## Project Bulk Modules to Single-cell Data

- Replace the discrete gene module assignments from sc.modules (UMAP+Louvain) to bulk.modules (WGNCA) from [here](https://rajlabmssm.github.io/WGCNA_MyND/2nd_batch_Analysis/protein_coding/04_sva/deepSplit_4/signed/04_page_sva_deep4_signed_PC.html#boxplot_mergedmes).  
- The "NA" module simply represents a collection of genes that were present in the sc dataset but not in the bulk dataset.

### PD vs. Control DGE Modules 

- Plot only the module that were differentially expressed (Wilcoxon nominal P-value < 0.05) between PD and controls in the bulk data.
- Also plot **darkmagenta**, which is of interest due to its enrichemnt in mitochondrial GO terms.

**Results Summary**  
- Bulk.modules differentially expressed between PD-Controls   
  + *Preferentially expressed in canonical monocytes*
    + darkmagenta (particularly the cells furthest from intermediate monocytes, suggesting these cells are the least activated)
    + pink  
    + honeydew1  
  + *Preferentially expressed in intermediate monocytes*  
    + purple (not clear separation but trending)


```r
merged.mods <- data.table:::merge.data.table(bulk.modules %>% dplyr::select(bulk.module=module, gene.symbol), 
                                             sc.modules, 
                                             by= "gene.symbol",
                                             all.y = T) %>%
  dplyr::select(gene.symbol, module=bulk.module, supermodule, dim_1, dim_2, dim_3) %>% unique() 



# Selected
select.modules <- c("lavenderblush3","darkturquoise","midnightblue","salmon4",
                    "honeydew1","pink","darkmagenta","green","blue","brown",
                    "darkolivegreen","darkorange2","purple","skyblue3",
                    "thistle1","thistle2") # sort(unique(merged.mods$module))[1:15]
p <- monocle3::plot_cells(cds,
                     genes=subset(merged.mods, module %in% select.modules),
                     group_cells_by="cluster",
                     color_cells_by="cluster",
                     show_trajectory_graph=F)
print(p)
```

<img src="scRNAseq_Monocle3_Modules_files/figure-html/Project Bulk Modules to Single-cell Data - PD.vs.Ctrl-1.png" width="960" />

```r
ggsave("./Results/sc_to_bulk.png", plot=p, dpi=400, width=10, height=7) 
```

### All Modules  

I also plotted all 60+ of the bulk.modules (ordered alphabetically) to see if any of them were preferentially expressed in a a given monocyte subtype. Darkmagenta remains one of the most striking examples, but there are other as well.

- Preferentially expressed in canonical monocytes  
  + grey80
  + ivory
- Preferentially expressed in intermediate monocytes  
  + firebrick4
  + orangered4
  + navajowhite2
  + red (trending)


```r
p2 <- monocle3::plot_cells(cds,
                     genes=merged.mods,
                     group_cells_by="cluster",
                     color_cells_by="cluster",
                     show_trajectory_graph=F)
print(p2)
```

<img src="scRNAseq_Monocle3_Modules_files/figure-html/Project Bulk Modules to Single-cell Data - All-1.png" width="1920" />

```r
ggsave("./Results/sc_to_bulk.all.png", plot=p2, dpi=400, width=20, height=14)
```

## Module-Module Overlap

- Calculate enrichment scores for all combinations of bulk RNAseq and scRNAseq gene modules.  
- This tells you whether some modules are detectable in both datasets.  


```r
genomeSize <- length(union(unique(bulk.modules$gene.symbol), monocle3::fData(cds)$gene_short_name))

MOD.dt <- module_vs_module_enrichment(mod.df1 = bulk.modules, 
                                      mod.df2 = sc.modules, 
                                      genomeSize = genomeSize, 
                                      save_file="./Results/bulk.modules_vs_sc.modules.txt",
                                      verbose = F) 
```

```
## [1] "mod.df1 contains 66 unique modules."
## [1] "mod.df2 contains 15 unique modules."
## [1] "Conducting enrichment tests on 990 module-module combinations..."
## [1] "990 enrichment tests conducted in 1.42 seconds."
## [1] "21 enrichment tests were significant (at FDR ≤ 0.05)."
## [1] "21.21% of mod.df1 modules showed enrichment for a mod.df2 module."
## [1] "86.67% of mod.df2 modules showed enrichment for a mod.df1 module."
```

```r
createDT(MOD.dt)  
```

preservef70e63eb39091065

### Plot Module-Module Overlap


```r
# library(heatmaply)
mod.cast <- reshape2::acast(MOD.dt, 
                            Module1.name~Module2.name, 
                            value.var="Jaccard.index")
hm <- heatmaply::heatmaply(mod.cast, 
                            k_row = 5, k_col = 5, 
                            height = 10, 
                            column_text_angle = 0, 
                            xlab = "scRNA-seq Louvain Modules",
                            ylab = "bulk RNA-seq WGCNA Modules", 
                            cexRow = .7,
                            key.title = "Jaccard Index") # file = "module-module.heatmap.html"
```

```
## Registered S3 method overwritten by 'seriation':
##   method         from 
##   reorder.hclust gclus
```

```r
# shiny::includeHTML("module-module.heatmap.html")
htmltools::tagList(setNames(list(hm), NULL)) 
```

preserve48b343353ca6ce4c


## Module GO Enrichment

- Iteratively test for enrichment of ontological terms in all modules (both bulk and single-cell).  
- Take a while for long lists of modules (e.g. the 60+ bulk-RNA-seq modules).


```r
mod.enrich.path <- "./Results/gprofiler2.module.enrichment.txt"
if(file.exists(mod.enrich.path)){
  print(paste("Module GO Enrichment file detected. Importing...",mod.enrich.path))
  mod.enrich <- data.table::fread(mod.enrich.path)
} else {
  bulk.mod.enrich <- gprofiler2.module_enrichment(bulk.modules, dataset = "bulk-RNA-seq")
  sc.mod.enrich <- gprofiler2.module_enrichment(sc.modules, dataset = "sc-RNA-seq")
  # Merge
  mod.enrich <- data.table::rbindlist(list(bulk.mod.enrich, sc.mod.enrich))
  # Correct for multiple tests
  mod.enrich <- mod.enrich %>% dplyr::mutate(FDR = p.adjust(p_value, method="fdr"), 
                                             Bonferroni = p.adjust(p_value, method="bonferroni")) %>% data.frame()
  mod.enrich[mod.enrich == "NULL"] <- NA # is.na(mod.enrich) <- mod.enrich == "NULL" 
  data.table::fwrite(mod.enrich, mod.enrich.path, sep="\t") 
}
```

```
## [1] "Module GO Enrichment file detected. Importing... ./Results/gprofiler2.module.enrichment.txt"
```

```r
mod.enrich.sig <- subset(mod.enrich, Bonferroni <= 0.05)
createDT(head(mod.enrich.sig)) 
```

preserve6cf05d91c0a61264

## Module-Module Similarity: Gene-level vs. Term-level 

- Compute how similar each of the top conserved module-module pairs are in terms of GO enrichment.  


```r
mm.overlap <- module_vs_module_enrichment.terms(MOD.dt = MOD.dt, 
                                                mod.enrich = mod.enrich.sig,  
                                                top.terms = F,
                                                row.limit = F)
```

```
## [1] "Comparing overlap of enrichment terms for all module-module combinations."
```

```r
createDT(mm.overlap)
```

preservefb2e709345cac9e3

### Histogram


```r
# Prepare dataframe
jaccard.corr <- merge(x = MOD.dt[,c("Module1.name","Module2.name","Jaccard.index")]  %>%
                        dplyr::rename(Jaccard.index.gene=Jaccard.index), 
                      y = mm.overlap[,c("Module1.name","Module2.name","Jaccard.index", 
                                        "Module1.top.term","Module2.top.term")] %>% 
                        dplyr::rename(Jaccard.index.term=Jaccard.index), 
                      by=c("Module1.name","Module2.name"))
# Remove zeros
# jaccard.corr <- subset(jaccard.corr, Jaccard.index.gene > 0 & Jaccard.index.term > 0)
## Histogram
ggplot(jaccard.corr) + 
  geom_histogram(aes(Jaccard.index.gene, fill="Jaccard.index.gene"), alpha=.5) + 
  geom_histogram(aes(Jaccard.index.term, fill="Jaccard.index.term"), alpha=.5)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<img src="scRNAseq_Monocle3_Modules_files/figure-html/Gene-level vs. Term-level - histogram-1.png" width="672" />

### Scatterplot

- See how well module-module similarity scores (Jaccard index) correlate between the gene-level comparisons vs. the term-level comparisons.  
- Varying how many of the top enrichment terms to use to calculate term-level module-module similarity (e.g. 10, 50, 100, all) did not seem to signficantly impact the correlation between methods (generaly hovered between Pearson's r = .4-.5).
<br>
- Label each pair of modules in the scatterplot using the top ontological enrichment terms for each module respectively. 
- Jaccard.sum is the sum of both the gene-level and term-level Jaccard similarity scores.  
  + This can be used as a way to identify module-module pairs that were similar at multiple levels (and thus have strong evidence of conservation).  


```r
## Correlation scatterplot
jaccard.corr <- jaccard.corr %>%
# label.dat <- subset(jaccard.corr, (Jaccard.index.gene >= .05| Jaccard.index.term >=.2) ) %>% 
  dplyr::mutate(Label=paste0("x: ",round(Jaccard.index.gene,2),
                             "\ny: ",round(Jaccard.index.term,2),
                             "\n\n",
                             # Mod 1
                             "Mod1 name: ",Module1.name,
                             "\nTop Term: ",Module1.top.term,
                             "\n\n", 
                             # Mod 2
                             "Mod2 name: ",Module2.name,
                             "\nTop Term: ",Module2.top.term),
                Jaccard.sum = Jaccard.index.gene + Jaccard.index.term)

corr.stat <- cor(jaccard.corr$Jaccard.index.gene, jaccard.corr$Jaccard.index.term, method="pearson")
# Scatterplot
sp <- ggplot(jaccard.corr, aes(x=Jaccard.index.gene, y=Jaccard.index.term, text=Label, color=Jaccard.sum)) + 
  geom_point() +
  # ggrepel::geom_text_repel(data = label.dat, aes(label=Label)) +
  geom_smooth(method=lm, se=T) + 
  labs(title = "Correlation between gene-level vs.term-level module-module Jaccard similarities",
       subtitle = paste("Pearson's r =",round(corr.stat,3)))
ggp <- plotly::ggplotly(sp, tooltip ="text" )
htmltools::tagList(list(ggp)) 
```

preserve320781d702203281


# Plot Bulk Hubs on Single-Cell Data


```r
# module-module comparisons
MOD.dt <- data.table::fread("./Results/bulk.modules_vs_sc.modules.txt") %>% 
  arrange(FDR)
top.bulk.mods <- head(MOD.dt, 10)$Module1.name %>% unique()
# One hub per module
bulk.hubs <- data.table::fread("./Data/Katia_Modules/hubs_gencode30.tbl.txt")
top.bulk.hubs <- subset(bulk.hubs, module %in% top.bulk.mods) 
# Plot hub expression in sc data
monocle3::plot_cells(cds,
                     genes=top.bulk.hubs$symbol,
                     group_cells_by="cluster",
                     color_cells_by="cluster",
                     show_trajectory_graph=F) #  rasterize = T (for super high-res figures) 
```

<img src="scRNAseq_Monocle3_Modules_files/figure-html/Plot Bulk Hubs on Single-Cell Data-1.png" width="672" />

 
 
 
