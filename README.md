# PD_scRNAseq 
Analysis pipelines for Parkinson's Disease scRNA-seq data from CD14+ peripheral monocytes.  

---  

## Seurat Results

The key difference with this set of analyses compared to Monocle2 below is that Seurat only tests a small subset a genes that meet certain criterion (e.g. are expressed in a minimum proportion of cells from each comparison group) to reduce the multiple-testing burden. This also makes sure that the DEGs that you do get are visibly different when plotting their expression in UMAP/PCA/t-SNE space.

### [Dimensionality Reduction + Clustering + Differential Gene Expression](https://rajlabmssm.github.io/PD_scRNAseq/scRNAseq_Seurat.html)

<hr><br>

## Monocle3 Results

### [Dimensionality Reduction + Pseudotime + Clustering + Cell Type Identification](https://rajlabmssm.github.io/PD_scRNAseq/scRNAseq_Monocle3.html)  

### [Differential Gene Expression + Top Specific Markers](https://rajlabmssm.github.io/PD_scRNAseq/scRNAseq_Monocle3_DGE.html)  
- Only top 1000 DGEs displayed in tables (page gets too slow otherwise). See **Full DGE Results Tables** link below for full results.  
- [Full DGE Results Tables](https://github.com/RajLabMSSM/PD_scRNAseq/tree/master/Results)  
 
### [DGE List Enrichment](https://rajlabmssm.github.io/PD_scRNAseq/scRNAseq_Monocle3_Enrich.html)  

### [Bulk vs. Single-cell Co-expression Modules](https://rajlabmssm.github.io/PD_scRNAseq/scRNAseq_Monocle3_Modules.html)  


---
### Created by:  
Brian M. Schilder, Bioinformatician II  
Raj Lab, Department of Neuroscience  
Icahn School of Medicine at Mount Sinai 
Nnew York City, New York
