import dataclasses
import json

from cogent3.app.composable import define_app
from cogent3.app.result import model_result
from cogent3.app.typing import SerialisableType
from cogent3.core.tree import PhyloNode
from cogent3.evolve.parameter_controller import AlignmentLikelihoodFunction

from phylo_limits.check_ident import validate_nodes
from phylo_limits.diagnose import diagonse
from phylo_limits.p_classifier import check_all_psubs
from phylo_limits.project import project


@dataclasses.dataclass(slots=True)
class PhyloLimitRec:
    """the record of phylogenetic limits"""
    unique_id: str 
    model_name: str
    tree: PhyloNode
    lf: AlignmentLikelihoodFunction
    psubs_class: dict 
    boundary_values: dict 
    bad_nodes: set
    identifiable: bool
    projection: list
        
    def to_rich_dict(self) -> dict:
        return {"id":self.unique_id, 
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
    bad_nodes = validate_nodes(lf=lf, strictly=strictly)
    identifiable= not bool(bad_nodes) # if there are any bad nodes show up, unidentifiable
    projection = project(lf) if identifiable else [] # only project rate matrix when the model is identifiable

    rec = PhyloLimitRec(unique_id=model_res.source, 
                        model_name=model_res.name,
                        tree=model_res.tree,
                        lf=lf,
                        psubs_class=check_all_psubs(lf=lf, strictly=strictly),
                        boundary_values=diagonse(lf=lf),
                        bad_nodes=bad_nodes,
                        identifiable= identifiable,
                        projection=projection)
    
    return json.dumps(rec.to_rich_dict())
