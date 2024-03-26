import py3Dmol
import stk, stko
from rdkit import Chem
from rdkit.Chem import AllChem
from ipywidgets import interact, interactive, fixed


def draw_mol(p,style='stick'):
        p.setStyle({style:{'colorscheme':'grayCarbon', 'singleBonds':True}})
        return p.show()

def show_molecule_from_file(filename):
    with open(filename,'r') as st:
       data = "".join(st.readlines())
    type = filename.split(".")[-1]
    p = py3Dmol.view(width=400,height=400)
    p.addModel(data,type)
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    styles = ["stick","sphere"]
    if type == "pdb":
        styles.append("cartoon")
    return interact(draw_mol, p=fixed(p), style=styles)

def show_molecule(stk_mol):
    data = AllChem.MolToMolBlock(stk_mol.to_rdkit_mol())
    p = py3Dmol.view(width=400,height=400)
    p.addModel(data,"mol")
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    styles = ["stick","sphere"]
    return interact(draw_mol, p=fixed(p), style=styles)