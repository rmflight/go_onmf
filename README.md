---
title: "RMF Investigation of GO-ONMF"
author: "Robert M Flight <rflight79@gmail.com>"
date: "2015-11-30 15:51:30"
output: md_document
---

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
double checked in `R` (see ONMF_source/brca/exp_onmf_brca.m)


```
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
```

## Prerequisites

To continue this analysis, you will need some other `R` packages installed.


```r
install.packages(c("readr", "dplyr"))
library(BiocInstaller)
biocLite(c("GO.db", "graph", "org.Hs.eg.db"))
```


```r
library(readr)
library(dplyr)
library(GO.db)
library(org.Hs.eg.db)
```

## Check Actual Input Data

Lets check that the actual data being used as input looks right.


```r
go_map <- read_csv("ONMF_source/brca/go_(merged).csv", col_names = c("loc", "GO"))
gene_map <- read_csv("ONMF_source/brca/gene_(merged).csv", col_names = c("loc", "GENE"))
gene2go <- read_csv("ONMF_source/brca/network_gene2go(merged).csv",
                    col_names = c("gene", "score", "GO"))
go2go <- read_csv("ONMF_source/brca/network_allgo2allgo(merged).csv",
                  col_names = c("GO1", "distance", "GO2"))
```

And translate to something we can actually use.


```r
gene2go$geneid <- unlist(gene_map[match(gene2go$gene, gene_map$loc), "GENE"],
                         use.names = FALSE)
gene2go$goid <- unlist(go_map[match(gene2go$GO, go_map$loc), "GO"],
                       use.names = FALSE)
```

And then query `org.Hs.eg.db` for the annotations as well.


```r
hs_gene2go <- select(org.Hs.eg.db, keys = unique(gene2go$geneid),
                     columns = c("GO", "ONTOLOGY"), keytype = "SYMBOL")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
hs_gene2go <- dplyr::filter(hs_gene2go, ONTOLOGY == "BP", EVIDENCE != "IEA")
```

And compare the annotations for the genes between them.


```r
annotation_overlap <- vapply(unique(hs_gene2go$SYMBOL), function(in_loc){
  hs_data <- dplyr::filter(hs_gene2go, SYMBOL == in_loc) %>%
    dplyr::select(., GO) %>% unlist()
  onmf_data <- dplyr::filter(gene2go, geneid == in_loc) %>%
    dplyr::select(., goid) %>% unlist()
  sum(hs_data %in% onmf_data) / length(hs_data)
}, numeric(1))
```

And check the distribution of overlap:


```r
hist(annotation_overlap)
```

![plot of chunk distribution_overlap](figure/distribution_overlap-1.png) 

That is great, over 4000 of the genes have what appear to be very similar terms!

Lets also check the **GO2GO** data.


```r
go_notself <- dplyr::filter(go2go, distance != 1)
go_notself$ID1 <- unlist(go_map[match(go_notself$GO1, go_map$loc), "GO"], use.names = FALSE)
go_notself$ID2 <- unlist(go_map[match(go_notself$GO2, go_map$loc), "GO"], use.names = FALSE)

onmf_pc <- split(go_notself$ID2, go_notself$ID1)
go_pc <- mget(names(onmf_pc), GOBPCHILDREN, ifnotfound = NA)

go2go_overlap <- vapply(names(onmf_pc), function(in_loc){
  go_data <- go_pc[[in_loc]]
  onmf_data <- onmf_pc[[in_loc]]
  sum(go_data %in% onmf_data) / length(go_data)
}, numeric(1))
```

And plot it:


```r
hist(go2go_overlap)
```

![plot of chunk plot_go2go_comparison](figure/plot_go2go_comparison-1.png) 



## Work With Single Gene

Now lets work with a single gene, and see what is going on. We'll keep it simple
and play with the first one in the list, *A2M*.




Now that we have the original and propagated indices in a format that we should
be able to read, we can compare them and check how things are being propagated,
whether to more or less specific terms.

### Read Indices


```r
org_indices <- scan(file = "ONMF_source/brca/org_notzero", what = numeric(), 
                    sep = " ")
org_indices <- org_indices[!is.na(org_indices)]
prop_indices <- scan(file = "ONMF_source/brca/prop_notzero", what = numeric(),
                     sep = " ")
prop_indices <- prop_indices[!is.na(prop_indices)]
```



### Actually Compare

First we need to actually get out the GO id's.


```r
org_go <- dplyr::filter(go_map, loc %in% org_indices) %>% dplyr::select(GO) %>% 
  unlist(., use.names = FALSE) %>% unique()
prop_go <- dplyr::filter(go_map, loc %in% prop_indices) %>% dplyr::select(GO) %>% 
  unlist(., use.names = FALSE) %>% unique()
```

And then what gene are we supposedly working with? I randomly chose the one with
the index of **783**.


```r
gene_id <- dplyr::filter(gene_map, loc %in% 783) %>% dplyr::select(GENE) %>% 
  unlist(., use.names = FALSE)
```

The gene name is CYP1B1.

Now, are the terms in `org_go` parents of the terms in `prop_go`? We will simply
make a calls to `GOBPANCESTORS` and save the results.


```r
# first check that everything is BP
org_terms <- GOTERM[org_go]
unique(Ontology(org_terms))
```

```
## [1] "BP"
```

```r
prop_terms <- GOTERM[prop_go]
unique(Ontology(prop_terms))
```

```
## [1] "BP"
```

```r
# get ancestors
prop_ancestors <- mget(prop_go, GOBPANCESTOR, ifnotfound = NA) %>% 
  unlist(., use.names = FALSE) %>% unique()

sum(org_go %in% prop_ancestors)
```

```
## [1] 8
```

Hmmm, there are 8 of the originals in the
propagated. This seems not right. Maybe they are children instead?


```r
org_ancestors <- mget(org_go, GOBPANCESTOR, ifnotfound = NA) %>%
  unlist(., use.names = FALSE) %>% unique()

sum(prop_go %in% org_ancestors)
```

```
## [1] 3
```

OK, this doesn't look right either. 

## Check GO Relationships

There is a file that has the distances between GO terms, encoded as:

* GO term 1 | distance | GO term 2

For adjacent terms, it appears that the distance is 0.5. So lets extract all
of them, and compare them to the data in `GO.db`.


```r
go2go <- read_csv("ONMF_source/brca/network_allgo2allgo(merged).csv",
                  col_names = c("GO1", "distance", "GO2"))
```




```r
go_child <- mget(names(onmf_parent_child), GOBPCHILDREN, ifnotfound = NA)
```

```
## Error in mget(names(onmf_parent_child), GOBPCHILDREN, ifnotfound = NA): error in evaluating the argument 'x' in selecting a method for function 'mget': Error: object 'onmf_parent_child' not found
```

And calculate the % overlap of the `go_child` to `onmf_parent_child`.


```r
perc_overlap <- vapply(names(onmf_parent_child), function(x){
  onmf_data <- unique(onmf_parent_child[[x]])
  go_data <- unique(go_child[[x]])
  
  sum(go_data %in% onmf_data) / length(go_data)
}, numeric(1))
```

```
## Error in vapply(names(onmf_parent_child), function(x) {: object 'onmf_parent_child' not found
```


```r
hist(perc_overlap)
```

```
## Error in hist(perc_overlap): object 'perc_overlap' not found
```

Some differences, but doesn't seem completely unexpected given that the data I'm
using is newer than the data they used.

## Check the GO Annotations

With this gene id, we can lookup the GO annotations for the gene stored in
`org.Hs.eg.db`.


```r
gene_go <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_id, columns = "GO",
                                 keytype = "SYMBOL")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
gene_go <- dplyr::filter(gene_go, ONTOLOGY == "BP", EVIDENCE != "IEA")
gene_allgo <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_id, columns = "GOALL",
                                    keytype = "SYMBOL")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
gene_allgo <- dplyr::filter(gene_allgo, ONTOLOGYALL == "BP", EVIDENCEALL != "IEA")
```

And lookup how many of the ONMF GO terms are in these two lists.


```r
sum(org_go %in% gene_go$GO) 
```

```
## [1] 1
```

