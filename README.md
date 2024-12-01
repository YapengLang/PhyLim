# phylim: a phylogenetic limit evaluation library built on [cogent3](https://cogent3.org/)

phylim evaluates the identifiability when estimating the phylogenetic tree using the Markov model. The identifiability is the key condition of the Markov model used in phylogenetics to fulfil consistency. Establishing identifiability relies on the organisation of five types of transition probability matrices on a phylogenetic tree. phylim provides a quick method to check the identifiability of a model fit, where we developed a main [cogent3 app](https://cogent3.org/doc/app/index.html), `check_identifiability`. 

The following content will demonstrate how to set up phylim and give some tutorials on `check_identifiability` and other associated apps.

## Installation

```pip install phylim```

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
>>> ident_check = get_app("phylim")
>>> record = ident_check(result)
>>> record.is_identifiable
True
```

Overall, `phylim` can wrap all such information up

```
>>> res.to_table()
```
|                           Phylo Limits Record                            |
|-------------------------------------------------------------------------|
| Source   | Model Name | Identifiable | Has Boundary Values | Version   |
|----------|------------|--------------|----------------------|-----------|
| Unknown  | GTR        | True         | True                 | 2024.9.20 |



<details>
<summary>Label all transition probability matrices in a model fit</summary>

You can call `classify_model_psubs` to give the category of all the matrices:

```
>>> from phylim import classify_model_psubs

>>> app = classify_model_psubs()
>>> labelled = app(result)
>>> labelled.to_rich_dict()
{'source': 'unknown', 'mcats': {(np.str_('Gorilla'),): 'DLC', (np.str_('Human'),): 'DLC', (np.str_('Mouse'),): 'DLC'}, 'version': '2024.9.20'}
```

</details>


<details>
<summary>Check if all parameter fits are within the boundary</summary>


```
>>> from phylim import check_fit_boundary

>>> app = check_fit_boundary()
>>> violations = app(result)
>>> violations.to_rich_dict()
{'source': 'unknown', 'vio': [{'par_name': 'C/T', 'init': np.float64(1.000000008361369e-06), 'lower': 1e-06, 'upper': 50}, {'par_name': 'A/T', 'init': np.float64(1.0000000181618708e-06), 'lower': 1e-06, 'upper': 50}], 'version': '2024.9.20'}
```

</details>
