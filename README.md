# GO ONMF

## Introduction

Kim and Yu, A Mutation Profile for Top-k Patient Search Exploiting Gene-Ontology
and Orthogonal Non-negative Matrix Factorization, Bioinformatics, 2015,
http://dx.doi.org/10.1093/bioinformatics/btv409

Kim and Yu recently published the above paper in which they develop a method
for creating a mutation-gene-ontology profile. As part of the overall procedure, 
they described (Section 2.3, *Gene Function Profile*):

> To reduce the term correlation, we use only the most specific terms, i.e. the leaf 
> node term after propagating the scores of the non-leaf terms down to the leaf node 
> terms. This approach also resolves the problem of evaluating genes annotated with 
> general term as the effect of the gene of function identification is spread out over 
> several leaf node terms.

This sounds troubling to my ears, because the Gene Ontology (GO) is a *directed 
acyclic graph* (DAG), wherein there are specific relationships between the terms, and
there is a directionality. There is also an expectation that when a gene product
is annotated to a *specific* term, then that gene product is also automatically
annotated to all of the *less-specific* **parent** terms in the DAG. 

What the procedure quoted above sounds like it is doing is essentially **adding**
more specific gene product to GO annotations, which as far as I know is not allowed.

This document is a record of my investigation as to whether *more specific* gene
product - GO annotations are generated during the creation of the mutation 
profiles.

## Software

Kim and Yu made an Octave version of their software available at https://sites.google.com/site/postechdm/research/implementation/orgos, and a
tarball of the Octave version can be downloaded from https://sites.google.com/site/postechdm/research/implementation/orgos/ONMF_octave-simple.zip

```
wget https://sites.google.com/site/postechdm/research/implementation/orgos/ONMF_octave-simple.zip
mkdir ONMF_source
unzip ONMF_octave-simple.zip -d ONMF_source
```

As part of the software package, a few data files are provided. The ones that
are important for our purposes include: 

* ONMF_source/brca/go_(merged).csv - the index to GO term file
* ONMF_source/brca/gene_(merged).csv - the index to gene file
* exp_onmf_brca.m - the file that runs the *brca* analysis
* network_gene2go(merged).csv - the gene 2 GO annotations

## Running Code

After modifying *exp_onmf_brca.m* with:

> library_path = '/home/rmflight/Projects/personal/onmf/ONMF_source/ONMF_octave';

And in `Octave` I ran:

```
cd brca
run exp_onmf_brca.m
```

As this runs, it generates *gene2go.mat*, which is the propogated scores of gene
to GO associations.

We will use `Octave` to generate some indices and print them to files that can be
double checked in `R`.

```
# from run_onmf_brca.m
gene2go = load('network_gene2go(merged).csv');
go2go = load('network_allgo2allgo(merged).csv');
patient2gene = getSparse(patient2gene);
gene2go = getSparse(gene2go);
go2go = getSparse(go2go);

patient2gene(find(patient2gene)) = 1;
gene2go(find(gene2go)) = 1;

# my own to save the data out
org_notzero = find(gene2go(783, :));
save -ascii org_notzero org_notzero;

# from ONMF_octave/TuneData.m
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

# save out the indices
prop_notzero = find(gene2go(783, :));
save -ascii prop_notzero prop_notzero;
```
