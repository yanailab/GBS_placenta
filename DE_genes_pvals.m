function [clust_genes,gene_pvals,clust_genes_names] = DE_genes(exp_mat,clust_vec,gene_names,p_thresh)
% perform ranksum test between each group vs all others
diff_mat = NaN(length(gene_names),max(clust_vec));
for f = 1:max(clust_vec)    
    for t = 1:length(gene_names)
        diff_mat(t,f) = ranksum(exp_mat(t,clust_vec == f),exp_mat(t,~(clust_vec == f)));
    end
end
significant_genes = diff_mat<p_thresh;
% calculate the sign of each gene between groups
sign_mat = NaN(length(gene_names),max(clust_vec));
for g = 1:max(clust_vec)
    for v = 1:length(gene_names)
    sign_mat(v,g) = mean(exp_mat(v,clust_vec==g)-mean(exp_mat(v,~(clust_vec==g)),2))>0;
    end
end
%select the top 100 genes for each (this could be a problem...depends how
%many u have)
cluster_genes = significant_genes&sign_mat;
[q,i] = sort(diff_mat + ~cluster_genes);

% top clust genes is the matrix of indices
clust_genes = NaN(100,max(clust_vec));
for m = 1:max(clust_vec)
    clust_genes(:,m) = i(1:100,m);
end

% give actual pvals
gene_pvals = NaN(100,max(clust_vec));
for g = 1:max(clust_vec)
    gene_pvals(:,g) = q(1:100,g);
end

%gives names in a mat
clust_genes_names = gene_names(clust_genes);
end
