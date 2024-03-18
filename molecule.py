from rdkit.Chem import AllChem
import py3Dmol

# Imports for visualization script
import stk, stko
from rdkit import Chem 
from rdkit.Chem import AllChem as rdkit
from collections import defaultdict
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import py3Dmol
from IPython.display import Image
import matplotlib.pyplot as plt
import subprocess
import time
import stk
import stko
import spindry as spd

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
