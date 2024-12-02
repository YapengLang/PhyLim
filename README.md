# phylim: a phylogenetic limit evaluation library built on [cogent3](https://cogent3.org/)

phylim evaluates the identifiability when estimating the phylogenetic tree using the Markov model. The identifiability is the key condition of the Markov model used in phylogenetics to fulfil consistency. Establishing identifiability relies on the organisation of five types of transition probability matrices on a phylogenetic tree. phylim provides a quick, handy method to check the identifiability of a model fit, where we developed a main [cogent3 app](https://cogent3.org/doc/app/index.html), `phylim`. phylim is compatible with [piqtree2](https://github.com/iqtree/piqtree2), a python library that exposes features from iqtree2.

The following content will demonstrate how to set up phylim and give some tutorials on `phylim` and other associated apps.

## Installation

```pip install phylim```

Let's see if it has been done successfully. In the package directory:

```pytest```

## Run the check of identifiability

If you fit a model to an alignment and get the model result:

```python
>>> from cogent3 import get_app, make_aligned_seqs

>>> algn = make_aligned_seqs(
...    {
...        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
...        "Gorilla": "ATGCGGCGCGCGGAGGCCGCGCTCGCGGAG",
...        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
...    },
...    info={"moltype": "dna", "source": "foo"},
... )

>>> app_fit = get_app("model", "GTR")
>>> result = app_fit(algn)
```

You can easily check the identifiability by:

```python
>>> app_ident_check = get_app("phylim")

>>> record = app_ident_check(result)
>>> record.is_identifiable
True
```

The `phylim` app wraps all information about phylogenetic limits.

```python
>>> record
```
</caption>
<thead class="head_cell">
<th>source</th><th>model name</th><th>identifiable</th><th>has boundary values</th><th>version</th>
</thead>
<tbody>
<tr><td><span class="c3col_left">brca1.fasta</span></td><td><span class="c3col_left">GTR</span></td><td><span class="c3col_left">True</span></td><td><span class="c3col_left">False</span></td><td><span class="c3col_left">2024.9.20</span></td></tr>
</tbody>
</table>

</div>

You can also use features like classifying all matrices or checking boundary values in a model fit.

<details>
<summary>Label all transition probability matrices in a model fit</summary>

You can call `classify_model_psubs` to give the category of all the matrices:

```python
>>> from phylim import classify_model_psubs

>>> labelled = classify_model_psubs(result)
>>> labelled.to_rich_dict()
{'source': 'foo', 'mcats': {(np.str_('Gorilla'),): 'DLC', (np.str_('Human'),): 'DLC', (np.str_('Mouse'),): 'DLC'}, 'version': '2024.9.20'}
```

</details>


<details>
<summary>Check if all parameter fits are within the boundary</summary>


```
>>> from phylim import check_fit_boundary

>>> violations = check_fit_boundary(result)
>>> violations.to_rich_dict()
{'source': 'foo', 'vio': [{'par_name': 'C/T', 'init': np.float64(1.000000008361369e-06), 'lower': 1e-06, 'upper': 50}, {'par_name': 'A/T', 'init': np.float64(1.0000000181618708e-06), 'lower': 1e-06, 'upper': 50}], 'version': '2024.9.20'}
```

</details>


## Colour the edges for a phylogenetic tree based on matrix categories


If you obtain a model fit, phylim can visualise the tree with labeled matrices. 

phylim provides an app, `phylim_colour_edges`, which takes an edge-matrix category map and colours the edges:

```python
>>> from phylim import classify_model_psubs

>>> edge_to_cat = classify_model_psubs(result)
>>> tree = result.tree
>>> app_colour_edge = get_app("phylim_colour_edges", edge_to_cat)
>>> app_colour_edge(tree)
```

![tree1](https://figshare.com/articles/figure/tree2_png/27936765)

You can also color edges with a user-defined edge-matrix category map: 

```python
>>> from phylim.classify_matrix import SYMPATHETIC, DLC

>>> edge_to_cat = {"A":SYMPATHETIC, "B":SYMPATHETIC, "C":DLC}

>>> app_colour_edge = get_app("phylim_colour_edges", edge_to_cat)
>>> app_colour_edge(tree)
```
![tree2](https://figshare.com/articles/figure/tree1_png/27936768?file=50903022)


## Check identifiability for piqtree2

phylim provides an app, `phylim_tree_to_likelihoodfunction`, which allows you to build the likelihood function from a piqtree2 output tree.

```python
>>> from piqtree2 import build_tree

>>> tree = build_tree(algn, model="GTR")
>>> app_inverter = get_app("phylim_tree_to_likelihoodfunction")

>>> result = app_inverter(tree)
>>> app_ident_check(result)
```
