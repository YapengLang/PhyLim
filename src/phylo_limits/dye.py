from cogent3.core.tree import PhyloNode
from cogent3.draw.dendrogram import Dendrogram
from cogent3.app.result import model_result
from phylo_limits.p_classifier import check_all_psubs

#TODO: colouring process lack exact groups
def colour_edges(model_res:model_result) -> Dendrogram:
    "a func accept a tree (in phynode ), output the colored edge based on the label dict"
    tree=model_res.tree
    psubs_dict=check_all_psubs(model_res)
    edge_Sympathetic=[]
    edge_I = []
    edge_D = []
    edge_C = []
    edge_L = []
    for i,j in psubs_dict.items():
        if j["class"] == "Sympathetic":
            edge_Sympathetic.append(i[0])
        elif j["class"] == "Identity":
            edge_I.append(i[0])
        elif j["class"] =="DLC":
            edge_D.append(i[0])    
        elif j["class"] == "Chainsaw":
            edge_C.append(i[0])
        elif j["class"] == "Limit":
            edge_L.append(i[0])
    
    fig = tree.get_figure()  
    fig.style_edges(edges = edge_Sympathetic, line=dict(color="blue", width=1), legendgroup="Sym")
    fig.style_edges(edges = edge_C, line=dict(color="red", width=1), legendgroup="Chainsaw")
    fig.style_edges(edges = edge_D, line=dict(color="grey", width=1), legendgroup="DLC")
    return fig