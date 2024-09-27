# phylo-limits: a phylogenetic limit evaluation library built on [cogent3](https://cogent3.org/)

phylo-limits evaluates the identifiability when estimating the phylogenetic tree using the Markov model. The identifiability is the key condition of the Markov model used in phylogenetics to fulfil consistency. Establishing identifiability relies on the organisation of five types of transition probability matrices on a phylogenetic tree. phylo-limits provides a quick method to check the identifiability of a model fit, where we developed a main [cogent3 app](https://cogent3.org/doc/app/index.html), `check_identifiability`. 

The following content will demonstrate how to set up phylo-limits and give some tutorials on `check_identifiability` and other associated apps.

## Installation

```pip install phylo-limits```

Let's see if it has been done successfully. In the package directory:

```pytest```

## Run the check of identifiability

If you fit a model to an alignment and get the model result:

```
>>> from cogent3 import get_app, make_aligned_seqs


>>> aln = make_aligned_seqs(
...    {
...        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
...        "Gorilla": "ATGCGGCGCGCGGAGGCCGCGCTCGCGGAG",
...        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
...    }
... )

>>> app_tr = get_app("model", "GTR")

>>> result = app_tr(aln)
```

Then, you can easily check the identifiability by:

```
>>> ident_check = get_app("check_identifiability")
>>> ident_check(result)
True
```

## Label all transition probability matrices in a model fit

You can call `classify_model_psubs` to give the category of all the matrices:

```
>>> label_psubs = get_app("classify_model_psubs")
>>> labelled = label_psubs(result)
>>> labelled.to_rich_dict()
{'source': 'unknown', 'mcats': {(np.str_('Gorilla'),): 'DLC', (np.str_('Human'),): 'DLC', (np.str_('Mouse'),): 'DLC'}, 'version': '2024.9.20'}
```

## Check if all parameter fits are within the boundary

```
>>> check_bound = get_app("check_fit_boundary")
>>> violations = check_bound(result)
>>> violations.to_rich_dict()
{'source': 'unknown', 'vio': [{'par_name': 'C/T', 'init': np.float64(1.000000008361369e-06), 'lower': 1e-06, 'upper': 50}, {'par_name': 'A/T', 'init': np.float64(1.0000000181618708e-06), 'lower': 1e-06, 'upper': 50}], 'version': '2024.9.20'}
```

## Overall, `generate_phylo_limit_record` can wrap all such information up

```
>>> check = get_app("generate_phylo_limit_record")
>>> record = check(result)
>>> record.to_rich_dict()
{'source': 'unknown', 'identifiability': True, 'strict': False, 'message': None, 'version': '2024.9.20', 'model_name': 'GTR', 'boundary_values': [{'par_name': 'C/T', 'init': np.float64(1.000000008361369e-06), 'lower': 1e-06, 'upper': 50}, {'par_name': 'A/T', 'init': np.float64(1.0000000181618708e-06), 'lower': 1e-06, 'upper': 50}], 'ISCL_mcats': {}}
```

