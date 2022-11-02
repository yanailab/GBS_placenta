function [info_genes,info_genes_names] = info_genes(exp_mat,max_thresh,fano_thresh,gene_names)
max_reads_per_cell = max(exp_mat')';
m =max_reads_per_cell>=max_thresh*mean(max_reads_per_cell);
fano_factor = var(exp_mat')'./mean(exp_mat,2);
f = fano_factor>fano_thresh*nanmean(fano_factor);
info_genes = f&m ;
info_genes_names = gene_names(info_genes);
end
