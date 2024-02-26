import json

from cogent3.app.composable import define_app
from cogent3.app.result import model_result
from cogent3.app.typing import SerialisableType

from phylo_limits.check_ident import validate_nodes
from phylo_limits.diagnose import diagonse
from phylo_limits.p_classifier import check_all_psubs
from phylo_limits.project import project


class PhyloLimitRec:
    """the record of phylogenetic limits"""
    def __init__(self, model_res:model_result):
        self.unique_id = model_res.source
        self.model_name = model_res.name
        self.tree = model_res.tree
        self.lf = model_res.lf

    def record(self, strictly=False):
        """record psubs classes, identifiability, boundary values and non-DLC projection 

        Args:
            "strictly" controls the sensitivity for Identity matrix (I); if false, treat I as DLC
        """    
        self.psubs_class = check_all_psubs(lf=self.lf, strictly=strictly)
        self.boundary_values = diagonse(self.lf)
        self.bad_nodes = validate_nodes(lf=self.lf, strictly=strictly)
        self.identifiable = not bool(self.bad_nodes) # if there are any bad nodes show up, unidentifiable
        self.projection = project(self.lf) if self.identifiable else [] # only project rate matrix when the model is identifiable
        
    def to_rich_dict(self) -> dict:
        return {"id":self.unique_id, 
                "model":self.model_name, 
                "tree":self.tree.get_newick(), 
                "psubs": self.psubs_class,
                "BVs":self.boundary_values,
                "identifiable":str(self.identifiable), 
                "bad_nodes":list(self.bad_nodes),
                "projection": self.projection}


@define_app
def app_for_all_infor(model_res:model_result) -> SerialisableType:
    rec = PhyloLimitRec(model_res=model_res)
    rec.record()
    return json.dumps(rec.to_rich_dict())
