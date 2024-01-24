import pytest
import numpy
from ylib.entro_meas import (
    multinomial_entropy,
    D_calculate_estimated_entropy,
    D_entropy_t_stat,
)

from cogent3 import load_aligned_seqs

from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage, importr
import rpy2.robjects as robjects

D_entropy_t_stat = D_entropy_t_stat()

adegenet = importr("adegenet")
obj = adegenet.fasta2DNAbin("../tests/data/primate_brca1.fasta")

obj2 = adegenet.fasta2DNAbin(
    "../tests/data/entropy_measure/100032_198959_515467_degap_nodege.fasta"
)

# move this into conftest.py
with open("../../entropy_saturation_test/saturation.index.R") as file:
    string = "".join(file.readlines())

stringr_c = SignatureTranslatedAnonymousPackage(string, "stringr_c")


@pytest.mark.parametrize(
    "n,p ",
    [
        (6, {"T": 0.2, "C": 0.1, "A": 0.3, "G": 0.4}),
        (10, {"T": 0.3, "C": 0.2, "A": 0.1, "G": 0.4}),
        (20, {"T": 0.25, "C": 0.25, "A": 0.3, "G": 0.2}),
    ],
)
def test_fss(n, p):
    result = multinomial_entropy(n, p)
    motif_probs = robjects.FloatVector(p.values())
    expected = stringr_c.multinomial_entropy(n=n, p=motif_probs)
    numpy.testing.assert_allclose(result, expected)


def test_exp():
    expected = numpy.asarray(stringr_c.calculate_index(obj)[0])[0]
    algn = load_aligned_seqs("../tests/data/primate_brca1.fasta", moltype="dna")
    result = D_calculate_estimated_entropy(algn)

    numpy.testing.assert_allclose(result, expected)


def test_tstat_case1():
    expected = numpy.asarray(stringr_c.calculate_index(obj)[3])[0]
    algn = load_aligned_seqs("../tests/data/primate_brca1.fasta", moltype="dna")
    result = D_entropy_t_stat(algn, only_varsites=False)[0]

    numpy.testing.assert_allclose(result, expected)


def test_tstat_case2():
    expected2 = numpy.asarray(stringr_c.calculate_index(obj2)[3])[0]
    algn2 = load_aligned_seqs(
        "../tests/data/entropy_measure/100032_198959_515467_degap_nodege.fasta",
        moltype="dna",
    )
    result2 = D_entropy_t_stat(algn2, only_varsites=False)[0]

    numpy.testing.assert_allclose(result2, expected2)
