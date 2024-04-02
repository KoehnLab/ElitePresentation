from ase.io import read, write
from ase.optimize import BFGS
from IPython.display import Markdown, display

class TrajectoryCalculator():
    
    def __init__(self, file, method = "GFN1-xTB", output = "trajectory.traj"):
        self.file = file
        self.method = method
        self.output = output
    
    def run_calculation(self):
        try:
            from tblite.ase import TBLite
        except:
            display(Markdown('<span style="color: #ff0000">Cannot run calculation. TBLite not provided.</span>.'))
            return
        atoms = read( self.file )
        atoms.calc = TBLite( method = self.method )
        dynamics = BFGS( atoms, trajectory = self.output )
        dynamics.run( fmax = 0.05 )