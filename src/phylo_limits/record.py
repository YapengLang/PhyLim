import dataclasses
import json

from cogent3.app.composable import define_app
from cogent3.app.result import model_result
from cogent3.app.typing import SerialisableType
from cogent3.core.tree import PhyloNode
from cogent3.evolve.parameter_controller import AlignmentLikelihoodFunction

from phylo_limits.check_boundary import diagonse
from phylo_limits.check_ident import has_valid_path
from phylo_limits.matrix_class import classify_psubs


@dataclasses.dataclass(slots=True)
class PhyloLimitRec:
    """the record of phylogenetic limits"""
    source: str 
    model_name: str
    tree: PhyloNode
    lf: AlignmentLikelihoodFunction 
    psubs_class: dict 
    boundary_values: dict 
    bad_nodes: set
    identifiable: bool
    projection: list
        
    def to_rich_dict(self) -> dict:
        return {"id":self.source, 
                "model":self.model_name, 
                "tree":self.tree.get_newick(), 
                "psubs_class": {str(k):v["class"] for k,v in self.psubs_class.items()},
                "BVs":self.boundary_values,
                "identifiable":str(self.identifiable), 
                "bad_nodes":list(self.bad_nodes),
                "projection": self.projection}


@define_app
def generate_record(model_res:model_result, strictly=False) -> SerialisableType:
    """record psubs classes, identifiability, boundary values and non-DLC projection 
    Args:
        "strictly" controls the sensitivity for Identity matrix (I); if false, treat I as DLC
    Return:
        string in json format 
    """  
    lf = model_res.lf
    bad_nodes = has_valid_path(lf=lf, strictly=strictly)
    identifiable= not bad_nodes # if there are any bad nodes show up, unidentifiable

    rec = PhyloLimitRec(source=model_res.source, 
                        model_name=model_res.name,
                        tree=model_res.tree,
                        lf=lf,
                        psubs_class=classify_psubs(lf=lf, strictly=strictly),
                        boundary_values=diagonse(lf=lf),
                        bad_nodes=bad_nodes,
                        identifiable= identifiable)
    
    return json.dumps(rec.to_rich_dict())
