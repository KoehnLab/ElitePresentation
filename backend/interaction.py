import py3Dmol
import stk, stko
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDetermineBonds
from ipywidgets import interact, fixed, widgets, GridspecLayout
from IPython.display import Markdown, display


def draw_mol(p,style='stick', pdb_color_scheme="nucleic"):
        if (style == "cartoon"):
            p.setStyle({style:{'colorscheme':pdb_color_scheme}})
        else:
            p.setStyle({style:{'colorscheme':'grayCarbon', 'singleBonds':True}})
        return p.show()

def show_molecule_from_xyz_file(filename):
    # Load molecule from XYZ file
    raw_mol = Chem.MolFromXYZFile(filename)
    mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineBonds(mol,charge=0)

    data = AllChem.MolToMolBlock(mol)
    p = py3Dmol.view(width=400,height=400)
    p.addModel(data,"mol")
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    styles = ["stick","sphere"]
    return interact(draw_mol, p=fixed(p), style=widgets.ToggleButtons(options=styles,disabled=False,button_style=''), pdb_color_scheme=fixed(""))


def show_molecule_from_file(filename):
    with open(filename,'r') as st:
       data = "".join(st.readlines())
    type = filename.split(".")[-1]
    p = py3Dmol.view(width=400,height=400)
    p.addModel(data,type)
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    styles = ["stick","sphere"]
    pdb_color_schemes = ["nucleic", "ssPyMol"]
    value = styles[0]
    if type == "pdb":
        styles.append("cartoon")
        value = styles[-1]
        return interact(draw_mol, p=fixed(p), style=widgets.ToggleButtons(options=styles,value=value,disabled=False,button_style=''), pdb_color_scheme=widgets.ToggleButtons(options=pdb_color_schemes,disabled=False,button_style='') )

    return interact(draw_mol, p=fixed(p), pdb_color_scheme=fixed(""), style=widgets.ToggleButtons(options=styles,value=value,disabled=False,button_style=''))

def show_molecule(stk_mol):
    data = AllChem.MolToMolBlock(stk_mol.to_rdkit_mol())
    p = py3Dmol.view(width=400,height=400)
    p.addModel(data,"mol")
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    styles = ["stick","sphere"]
    return interact(draw_mol, p=fixed(p), style=widgets.ToggleButtons(options=styles,disabled=False,button_style=''), pdb_color_scheme=fixed(""))

def __show_smiles(smiles = "CCCC", style = "stick", optimize = False):
    try:
        stk_mol = stk.BuildingBlock( smiles )
    except:
        display(Markdown('<span style="color: #ff0000">Ungültige SMILES Spezifikation</span>.'))
        return
    if optimize:
        stk_mol = stko.MMFF().optimize(stk_mol)
    data = AllChem.MolToMolBlock(stk_mol.to_rdkit_mol())
    p = py3Dmol.view(width=400,height=400)
    p.addModel(data,"mol")
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    return draw_mol(p, style)


def __show_structure_formula(stk_mol):
    formula = stk_mol.to_rdkit_mol()
    formula = Chem.RemoveHs( formula )
    AllChem.Compute2DCoords( formula )
    pil = Draw.MolToImage(formula, size=(400,400))
    return display(pil)

def __show_structure_formula(smiles = "CCCC", optimize = False):
    try:
        stk_mol = stk.BuildingBlock( smiles )
    except:
        display (Markdown('<span style="color: #ff0000">Ungültige SMILES Spezifikation</span>.'))
        return
    if optimize:
        stk_mol = stko.MMFF().optimize(stk_mol)
    formula = stk_mol.to_rdkit_mol()
    formula = Chem.RemoveHs( formula )
    AllChem.Compute2DCoords( formula ) 
    pil = Draw.MolToImage(formula, size=(400,400))
    return display(pil)


def __show_both(smiles = "CCCC", style = "stick", optimize = False):
    if isinstance(smiles, str):
        sm = widgets.Text(value = smiles, placeholder = 'String', description = 'SMILES', disabled = False, layout={'width': 'auto'})
    else:
        sm = widgets.Dropdown(options = smiles, value = smiles[0], description = "SMILES", layout={'width': 'auto'})
    
    st = widgets.ToggleButtons(
    options=style,
    #description='Render Styles:',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    )
    #st = widgets.RadioButtons(options = style, value = style[0], description = "Styles", layout={'width': 'max-content'})
    o = widgets.Checkbox(value = False, description = "Optimierung")

    out1 = widgets.interactive_output(__show_smiles, {'smiles': sm, 'style': st, 'optimize': o})
    out2 = widgets.interactive_output(__show_structure_formula, {'smiles': sm, 'optimize': o})

    grid = GridspecLayout(3, 4, height='400px')
    grid[0, :] = sm
    grid[1, 0] = st
    grid[1, 1] = o
    grid[2, 0] = out1
    grid[2, 1] = out2
    # display(HBox([sm,st,o]))
    # display(HBox([out1, out2]))
    display(grid)


def show_molecule_from_smiles(smiles = "CCCC", show_structure_formula = False):
    styles = ["stick","sphere"]
    if show_structure_formula:
        return __show_both(smiles = smiles, style = styles, optimize = [True, False])
    
    return interact(__show_smiles, smiles = smiles, style = widgets.ToggleButtons(options=styles,disabled=False,button_style=''), optimize = widgets.Checkbox(value = False, description = "Optimierung"))

def show_structure_from_file(file_path):
    mol = Chem.MolFromMolFile(file_path)
    mol = Chem.RemoveHs( mol )
    AllChem.Compute2DCoords( mol ) 
    pil = Chem.Draw.MolToImage(mol, size=(400,400))
    return display(pil)