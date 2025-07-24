"""A library for checking the limits of phylogenetic tree estimation."""

from phylim.apps import check_fit_boundary, classify_model_psubs
from phylim.classify_matrix import CHAINSAW, DLC, IDENTITY, LIMIT, SYMPATHETIC


check_fit_boundary = check_fit_boundary()
classify_model_psubs = classify_model_psubs()
