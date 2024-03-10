import py3Dmol
from ase.io.trajectory import Trajectory
from ase.io import read, write
from molecule import *

class TrajectoryVisualizer:
    
    view = py3Dmol.view(width=400, height=300)

    def __init__(self, file):
        self.file = file

    def create_trajectory(self):

        trajectory = Trajectory(self.file)

        molecules = []

        for atoms in trajectory:

            write( 'intermediary.pdb', atoms )

            with open( 'intermediary.pdb' ) as ifile: molecule = Molecule(ifile)

            molecules.append( molecule )

        for mol in molecules:
        
            for at in mol:
                
                if at["resname"] == "PRO":
                    
                    at["pymol"] = {"stick": {'color': "yellow"}}
                    
                elif at["resname"] == "GLY":
                    
                    at["pymol"] = {"stick": {'color': 'blue'}}

        models = ""

        for i, mol in enumerate(molecules):
            
            models += "MODEL " + str(i) + "\n"
            
            models += str(mol)
            
            models += "ENDMDL\n"
            
        self.view.addModelsAsFrames(models)

        return molecules

    def show_trajectory(self):

        molecules = self.create_trajectory()

        for i, at in enumerate(molecules[0]):
        
            default = {'stick':{'colorscheme':'greyCarbon'}}
            
            self.view.setStyle({'model': -1, 'serial': i+1}, at.get("pymol", default))

        self.view.zoomTo()

        self.view.animate({'loop': "forward"})

        self.view.show()