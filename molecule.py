from rdkit.Chem import AllChem
import py3Dmol
import stk, stko
from rdkit import Chem 
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

class Molecule:
    cyan_colorscheme = False
    def set_cyan_colorscheme(set):
        Molecule.cyan_colorscheme = set

def show_stk_mol(stk_mol):
    data = AllChem.MolToMolBlock(stk_mol.to_rdkit_mol())
    colorscheme = 'cyanCarbon' if Molecule.cyan_colorscheme else 'greyCarbon'
    
    p = py3Dmol.view(
        data=data,
        style={'stick':{'colorscheme':colorscheme, 'singleBonds':True}}, 
        width=400,
        height=400,
    )

    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    p.show()
