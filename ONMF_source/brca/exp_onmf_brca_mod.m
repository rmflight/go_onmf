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

org_notzero = find(gene2go(783, :));
save -ascii org_notzero org_notzero;

[baseSMData, patient2gene, gene2go_new, go2go_new] = TuneData(patient2gene, gene2go, go2go);

isequal(go2go_new, go2go)

prop_notzero = find(gene2go_new(783, :));
save -ascii prop_notzero prop_notzero;
