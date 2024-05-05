import dataclasses


@dataclasses.dataclass(slots=True)
class PhyloLimitRec:
    """the record of phylogenetic limits"""

    source: str
    model_name: str
    tree: str
    psubs_class: dict
    boundary_values: dict
    identifiable: bool
    bad_nodes: set

    def to_rich_dict(self) -> dict:
        return {
            "source": self.source,
            "model_name": self.model_name,
            "tree": self.tree,
            "psubs_class": self.psubs_class,
            "boundary_values": self.boundary_values,
            "identifiable": self.identifiable,
            "bad_nodes": self.bad_nodes,
        }
