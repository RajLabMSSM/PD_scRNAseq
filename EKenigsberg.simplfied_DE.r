# umitab - UMI table (genes, cells)
# mask_bg - A vector containing the cell_IDs of set1
# mask_fg - A vector containing the cell_IDs of set2
# nmin_umi_thresh - A threshold deterimning whether a gene will be considered expressed in a given cell for the gene filtering criterion
# nmin_cells_with_min_umi - A threshold deterimning the minimum number of expressing cells. Genes expressed in more cells will be included in the analusis.
# reg - regularization for the log fold change values
# nchunks - number of iterations
# n_per_chunk - number of permutations per chunk

DE_between_two_sets=function(umitab,mask_bg,mask_fg,nmin_umi_thresh=0,nmin_cells_with_min_umi=20,reg=1e-5,nchunks=100,n_per_chunk=1000){
  
  u=umitab[,c(mask_bg,mask_fg)]
  n_cells_with_numis_above_thresh=rowSums(u>nmin_umi_thresh)
  names(n_cells_with_numis_above_thresh)=rownames(u)
  gene_mask=n_cells_with_numis_above_thresh>nmin_cells_with_min_umi
  u=u[gene_mask,c(mask_bg,mask_fg)]
  message("Testing ",sum(gene_mask)," genes")
 
  obs_s=pmax(rowSums(u),0)
  obs_s_bg=pmax(rowSums(u[,mask_bg,drop=F]),0)
  obs_s_fg=pmax(obs_s-obs_s_bg,0)
  obs_m_bg=obs_s_bg/sum(obs_s_bg)
  obs_m_fg=obs_s_fg/sum(obs_s_fg)
  obs_log2_fc=log2((reg+obs_m_fg)/(reg+obs_m_bg))
  ncounts_bigger=rep(0,nrow(u))
  n1=length(mask_bg)
  mat=matrix(c(rep(T,n1),rep(F,ncol(u)-n1)),n_per_chunk,ncol(u),byrow =T)
  
  ntot=ncol(u)
  for (i in 1:nchunks){
    message(i)
    
    print(system.time({
      s_bg=pmax(u%*%apply(mat,1,sample,ntot),0)
    }))
    s_fg=pmax(obs_s-s_bg,0)
    m_bg=t(t(s_bg)/colSums(s_bg))
    m_fg=t(t(s_fg)/colSums(s_fg))
    log2_fc=log2((reg+m_fg)/(reg+m_bg))
    ncounts_bigger=ncounts_bigger+rowSums(abs(log2_fc)>=abs(obs_log2_fc))
    
    #    save(file="de_current_iter.rd",i)
  }
  p.value=ncounts_bigger/(i*n_per_chunk)
  adj.p.value=p.adjust(p.value,method = "BH")
  de_res=data.frame(counts_bg=obs_s_bg,counts_fg=obs_s_fg,freq_bg=obs_m_bg,freq_fg=obs_m_fg,n_cells_with_min_umis_above_thresh=n_cells_with_numis_above_thresh[gene_mask],log2_FC=obs_log2_fc,p.value=p.value,adj.p.value=adj.p.value)
  
  return(de_res)
}

