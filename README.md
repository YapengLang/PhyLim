# phylo_limits
This repo is a library for Yapeng Masters project scripts [see my thesis].

There is a main function -- `record_results.generate_record`, which works on a cogent3 model_result object. It will:

1. read the phylogeny and all p matrices
2. label all p matrices by classes 
3. check the identifiability based on theories 
4. check if there are boundary values in fits

