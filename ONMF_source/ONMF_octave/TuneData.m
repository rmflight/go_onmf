function [baseSMData, patient2gene, gene2go, go2go] = TuneData(patient2gene, gene2go, go2go, min_mutation)

[go_col] = size(gene2go, 2);
[go_row] = size(go2go, 1);

if (go_row ~= go_col)
    min_gonum = min(go_col, go_row);
    go2go = go2go(1:min_gonum, 1:min_gonum);
    gene2go = gene2go(:, 1:min_gonum);
end

[gene_col] = size(patient2gene, 2);
[gene_row] = size(gene2go, 1);
if (gene_col ~= gene_row)
    min_genenum = min(gene_row, gene_col);
    patient2gene = patient2gene(:, 1:min_genenum);
    gene2go = gene2go(1:min_genenum, :);
end

baseSMData.gene_indiv_mat = patient2gene;
if ~exist('gene2go.mat')
    gene2go = networkPropagateWithTof(gene2go, go2go, 1.0, 100);
    save('gene2go.mat', 'gene2go');
else
    load('gene2go.mat');
end
baseSMData.gene_indiv_mat = baseSMData.gene_indiv_mat * gene2go;
baseSMData.gene_indiv_mat(:, sum(baseSMData.gene_indiv_mat, 1) == 0) = [];
baseSMData.gene_indiv_mat = full(baseSMData.gene_indiv_mat);
% baseSMData.gene_indiv_mat = quantilenorm(baseSMData.gene_indiv_mat')';
baseSMData.gene_indiv_mat = diag(min(1./baseSMData.gene_indiv_mat')) * baseSMData.gene_indiv_mat;
