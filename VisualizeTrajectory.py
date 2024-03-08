import py3Dmol

class TrajectoryVisualizer:
    
    view = py3Dmol.view(width=400, height=300)

    def __init__(self, mols):
        self.molecules = mols

    def create_trajectory(self):

        for mol in self.molecules:
        
            for at in mol:
                
                if at["resname"] == "PRO":
                    
                    at["pymol"] = {"stick": {'color': "yellow"}}
                    
                elif at["resname"] == "GLY":
                    
                    at["pymol"] = {"stick": {'color': 'blue'}}

        models = ""

        for i, mol in enumerate(self.molecules):
            
            models += "MODEL " + str(i) + "\n"
            
            models += str(mol)
            
            models += "ENDMDL\n"
            
        self.view.addModelsAsFrames(models)

    def show_trajectory(self):

        self.create_trajectory()

        for i, at in enumerate(self.molecules[0]):
        
            default = {'stick':{'colorscheme':'greyCarbon'}}
            
            self.view.setStyle({'model': -1, 'serial': i+1}, at.get("pymol", default))

        self.view.zoomTo()

        self.view.animate({'loop': "forward"})

        self.view.show()