import json

from cogent3.app.composable import define_app
from cogent3.app.result import model_result
from cogent3.app.typing import SerialisableType

from phylo_limits.check_ident import check_ident
from phylo_limits.diagnose import diagonse
from phylo_limits.project import project


class PhyloLimitRec:
    """the record of phylogenetic limits"""
    def __init__(self, model_res:model_result):
        self.unique_id = model_res.source
        self.model_name = model_res.name
        self.tree = model_res.tree
        self.psub_dict = model_res.lf.get_all_psubs()
        self.qsub_dict = model_res.lf.get_all_rate_matrices(calibrated=False)
        self.qcsub_dict = model_res.lf.get_all_rate_matrices(calibrated=True)
        self.motif_probs = model_res.lf.get_motif_probs_by_node()
        self.lengths = model_res.lf.get_param_value

    def record(self):    
        self.identifiability = check_ident(...)
        self.boundary_values = diagonse(...)
        if self.identifiability:
            self.projection = project(...)
        else:
            self.projection = None 
        
    def to_rich_dict(self) -> dict:
        return print("")


@define_app
def app_for_all_infor(model_res:model_result) -> SerialisableType:
    rec = PhyloLimitRec(model_res=model_res)
    rec.record()
    return json.dumps(rec.to_rich_dict())
