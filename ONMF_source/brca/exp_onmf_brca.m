%% Orthogonal non-negative matrix factorization/Gene-Ontology-based stratification
library_path = '/home/rmflight/Projects/personal/onmf/ONMF_source/ONMF_octave';
addpath(genpath(library_path))

%% Convert somatic mutation cancer data to my data
patient2gene = load('network_patient2gene(merged).csv');
gene2go = load('network_gene2go(merged).csv');
go2go = load('network_allgo2allgo(merged).csv');
patient2gene = getSparse(patient2gene);
gene2go = getSparse(gene2go);
go2go = getSparse(go2go);

patient2gene(find(patient2gene)) = 1;
gene2go(find(gene2go)) = 1;

baseSMData = TuneData(patient2gene, gene2go, go2go);

fid = fopen('patient_patient2gene_row.csv');
str = textscan(fid, '%d %s', 'delimiter', ',');
baseSMData.sample_id = regexprep(str{2}, '(TCGA-..-....).*','$1');
fclose(fid);

%% Load phenotype
LoadBRCA;
[~,ia,ib] = intersect(baseSMData.sample_id, ppheno.sample_id);
clinical_data = [ppheno.days_survival(ib) ...
    ppheno.diag_age(ib) ...
    grp2idx(ppheno.estrogen_status(ib)) ...
    grp2idx(ppheno.progesterone_status(ib)) ...
    grp2idx(ppheno.hist_type(ib)) ...
    grp2idx(ppheno.margin_status(ib)) ...
    ppheno.survival(ib)];

%% Run NMF
ornmf_options = statset('nnmf');
ornmf_options.TolX = 1e-4;
ornmf_options.TolFun = 1e-4;
ornmf_options.MaxIter = 1500;

k = 4;
[W, Y, d] = ornmf(baseSMData.gene_indiv_mat', k, 'algorithm', 'mult', 'options', ornmf_options);
[~, subtype_label] = max(Y,[],1);

h1 = subtype_label(ia)';
h2 = clinical_data(:,3);
[table, l_x, l_y] = crosstable(h1, h2);
[p, chi2, df] = chisquare_test_independence(table);
fprintf('%d-ONMF: chi2:%f, p:%e\n',k,chi2,p);

clear baseSMData;

