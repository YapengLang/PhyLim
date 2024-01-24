# this file implements the entropy measure from Xia and Duchene.
# two functions only realise 1) the entropy on one site under full sat. 2) the observed entropy for given data
# however, the critical value for the following t-test is not provided in this file; this should be calced by duchene's
# and then turn back to do the t-test
from numpy import log, math, arange, array, partition
from scipy.special import comb
from scipy.stats import multinomial, ttest_1samp
from math import log as mathlog

from cogent3.core.alignment import ArrayAlignment
from cogent3.app.typing import AlignedSeqsType
from cogent3.app.composable import define_app
from cogent3 import get_moltype

fc = math.factorial


def multinomial_entropy(n_taxa: int, motif_probs: dict) -> float:
    """the entropy on one site under full sat, given number of taxa and motif probs"""
    x = list(arange(0, n_taxa))
    k = ["T", "C", "A", "G"]

    term3 = []
    for i in k:
        for j in x + [n_taxa]:
            each = (
                comb(n_taxa, j)
                * (motif_probs[i] ** j)
                * ((1 - motif_probs[i]) ** (n_taxa - j))
                * mathlog(math.factorial(j))
            )
            term3.append(each)

    return (
        -mathlog(fc(n_taxa))
        - n_taxa * sum(motif_probs[i] * mathlog(motif_probs[i]) for i in k)
        + sum(term3)
    )


def D_calculate_estimated_entropy(algn: ArrayAlignment, stepwise=False) -> float:
    """the observed entropy for given the whole alignment.
    here, duchene used the global motif probs as params the multinomial.
    so the information on one site was computed with assumption that the composition behind same to the global
    """
    DNA = get_moltype("dna")
    length = len(algn)
    n = algn.num_seqs
    motif_probs = algn.get_motif_probs()

    infors = []

    for i in range(length):
        counts = algn[i].counts().to_dict()
        prob = multinomial.pmf(
            x=[counts[b] for b in DNA],
            n=n,
            p=[motif_probs[b] for b in DNA],
        )
        infor = -log(prob)
        infors += [infor]
    return infors if stepwise else sum(infors) / length


def only_vars_indices(algn):
    counts_per_site = algn.counts_per_pos().to_dict()

    indices = []
    for pos, count in counts_per_site.items():
        count_list = array(list(count.values()))
        if count_list.max() < algn.num_seqs - 1 and partition(count_list, -2)[-2] > 1:
            indices.append(pos)
    return indices


@define_app
def D_entropy_t_stat(algn: ArrayAlignment, only_varsites=False) -> float:
    """construct t-statistic as per Duchene et al 2022"""
    algn_ndege = algn.no_degenerates(allow_gap=False)
    if only_varsites:
        if algn.num_seqs < 4:
            raise ValueError("Calc on infor sites only supports num_seqs > 3")
        algn_ndege = algn_ndege.take_positions(only_vars_indices(algn_ndege))
    stepwise_infor = D_calculate_estimated_entropy(algn=algn_ndege, stepwise=True)
    fss_entropy = multinomial_entropy(
        n_taxa=algn_ndege.num_seqs, motif_probs=algn_ndege.get_motif_probs()
    )
    return ttest_1samp(a=stepwise_infor, popmean=fss_entropy, alternative="two-sided")


# these used in cli for model fitting
@define_app
def only_varsites(val: AlignedSeqsType) -> AlignedSeqsType:
    """a composable app for obtaining infor sites only, as Duchene et al 2022 defined"""
    algn_ndege = val.no_degenerates(allow_gap=False)
    return algn_ndege.take_positions(only_vars_indices(algn_ndege))


@define_app
def de_gap(val: AlignedSeqsType) -> AlignedSeqsType:
    """degap"""
    return val.no_degenerates(allow_gap=False)
