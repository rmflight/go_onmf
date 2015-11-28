An Octave implementation of "A Mutation Profile for Top-k Patient Search Exploiting Gene-Ontology and Orthogonal Non-negative Matrix Factorization) (available in octave)

Authors: 
    Sungchul Kim (subright@postech.ac.kr)
    Lee Sael (sael@cs.stonybrook.edu)
    Hwanjo Yu (hwanjoyu@postech.ac.kr)

Abstract:
    Given a large quantity of genome mutation data collected from clinics, how can we search for similar patients? Similarity search based on patient mutation profiles can solve various translational bioinformatics tasks, including prognostics and treatment efficacy predictions for better clinical decision making through sheer volume of data. However, this is a challenging problem due to heterogeneous and sparse characteristics of the mutation data as well as its high dimensionality. To solve this problem we introduce a compact representation and search strategy based on Gene-Ontology and orthogonal non-negative matrix factorization. Results show that our method is able to identify and characterize clinically meaningful tumor subtypes better than the recently introduced Network Based Stratification method while enabling real-time search. To the best of our knowledge, this is the first attempt to simultaneously characterize and represent somatic mutational data for efficient search purposes.

Citation: 


--------------------------------------------------------------------------------------

- Tested platform:
	Linux & Mac running Octave	

- How to execute demo files for Stratifying mutation data (represented in Gene-Ontology) via Orthogonal Non-negative matrix factorization:
    * {cancer type} = {ucec, ov, luad, gbm, brca}
	1. Unpack project and data files. 
    2. Add library and data directories to your matlab path by changing the "library_path" in the demo files located in exp_onmf_{cancer type}.m files of each data directory: 
        > library_path = '../ONMF_octave'; 
        to 
        > library_path = 'location of "ONMF_octave" directory'; 
    3. Run Octave
        > octave
    4. Move to a {cancer type} directory
        octave:1> cd {cancer type}
    5. Run exp_onmf_{cancer type}.m
        octave:2> run exp_onmf_{cancer type}.m
       Outputs: "[ONMF] clusters: %d, chi-square:%f, p-value:%e" 

- To apply ONMF for your data, you need to prepare the three  data matrices and use the two functions described bellow:
    + Data Matrices
        1. patient2gene: a binary sparse matrix of patient-by-gene that indicates the occurrences of somatic mutation on genes for each patients
        2. gene2go: a binary sparse matrix of gene-by-go(or go-terms) that indicates the associtaion of gene or gene products with GO terms
        3. go2go: a sparse matrix of go-by-go that indicates the patient-child relationship between GO terms on Gene-Ontology.	

    + TuneData.m
        1. Syntax: [DataCell] = TuneData(patient2gene, gene2go, go2go);
	    2. Description: [DataCell] = TuneData(patient2gene, gene2go, go2go) returns the data structure by  somatic mutation data 
    + ornmf.m
        1. Syntax: [W, Y, d] = norm_ornnmf(X, k, 'algorithm', 'mult', 'options', options, ...);
	    2. Description: [W, Y, d] = norm_ornnmf is derived from 'nnmf.m'. It factors the non-negative n-by-m matrix X into non-negative factors W (n-by-k) and H (k-by-m). W consists of coefficient vectors to reconstruct the input matrix, and H consists of basis vectors.
	    3. Parameters: 'algorithm' should be 'mult' and options are derived from statset('nnmf'). 		
	
- Note
	This implementation omits pre-processing steps and ignores comparison with 'NBS' for simplification (i.e., mapping genes and patients between data used in NBS and ONMF), thus there are some variations on the result compared to the originial implementation.
