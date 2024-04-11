import torch
import schnetpack as spk
import schnetpack.transform as trn
import os
from ase.io import write
from ase.md.langevin import Langevin
from ase import Atoms, units
import py3Dmol
from IPython.display import Markdown, display
from ipywidgets import IntProgress
from rdkit import Chem
from ase.io.trajectory import Trajectory


class MLSimulation():

    def __init__(self, mol_file = "salicylic_acid.mol", model_file = "salicylic_acid_parameters"):
        self.mol_file = os.path.join( "data", mol_file )
        self.model_file = os.path.join( "data", model_file )
        self.device = torch.device( 'cuda' )
        self.script_dir = os.path.abspath(os.path.dirname(__file__))
        self.output_trajectory = "data/salicylic_acid_trajectory.traj"
        self.output_trajectory_xyz = "data/salicylic_acid_trajectory.xyz"
        self.sim_steps = 100
        self.sim_step_size = 10


    def run_simulation(self):
        # Read with rdkit
        mol = Chem.MolFromMolFile(self.mol_file, removeHs=False)
        # Convert for use in ase
        positions = mol.GetConformer().GetPositions()  # Get atomic positions
        atomic_numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]  # Get atomic numbers

        # Create ase atoms construct
        self.atoms = Atoms(numbers=atomic_numbers, positions=positions)

        # Simulation definitions
        calculator = spk.interfaces.SpkCalculator(
                model_file = self.model_file,
                neighbor_list = trn.ASENeighborList( cutoff = 5. ),
                energy_key = 'energy',
                force_key = 'forces',
                energy_unit = 'kcal/mol',
                position_unit = 'Ang',
                )
        self.atoms.calc = calculator


        dynamics = Langevin(
            self.atoms,
            timestep=1.0 * units.fs,
            temperature_K=298.0,  # temperature in K
            friction=0.01 / units.fs,
            trajectory = self.output_trajectory
        )

        # Init progress bar
        display(Markdown('Start der Simulation...'))
        f = IntProgress(min=0, max=self.sim_steps, layout={'width': 'auto'})
        display(f) # display the bar

        # Run the Simulation
        for i in range(self.sim_steps):
            dynamics.run( self.sim_step_size )
            f.value += 1


    def create_xyz_trajectory(self):
        traj = Trajectory( self.output_trajectory )

        path = "data/ml_scratch"
        if not os.path.exists(path): os.makedirs(path)
        
        # Add another progress bar
        display(Markdown('Konvertierung der Daten...'))
        f = IntProgress(min=0, max=len(traj), layout={'width': 'auto'})
        display(f)

        #write each structure from the .traj file in .xyz format
        for i in range(0, len(traj)):
            atoms = traj[i]
            string = 'structure%03d' % (i,) +'.xyz'
            outStruct = os.path.join(path, string)
            write(outStruct, atoms)
            #combines all optimization structures in one trajectory 
            #file
            inFile = open(os.path.join(path, 'structure%03d' % 
                        (i,)  +'.xyz'), 'r')
            fileStr = inFile.read()
            outFile = open(self.output_trajectory_xyz, 'a')
            outFile.write(fileStr)
            f.value += 1


    def show_result(self):
        if not os.path.exists(self.output_trajectory):
            display(Markdown('<span style="color: #ff0000">Keine Ergebisdatei gefunden. Bitte zuerst die Simulation laufen lassen.</span>.'))
            return
        
        self.create_xyz_trajectory()

        xyz_data = open(self.output_trajectory_xyz,"r").read()

        view = py3Dmol.view(width=400, height=300)
        view.addModelsAsFrames(xyz_data)
        view.animate({"loop": "forward"})
        # visualize with the stick option - can also consider spheres and more
        view.setStyle({"stick": {'colorscheme':'greyCarbon'}})
        view.show()


