from ase.io import read, write
from tblite.ase import TBLite
from ase.optimize import BFGS


class TrajectoryCalculator():
    
    def __init__(self, file, method = "GFN1-xTB", output = "trajectory.traj"):
        self.file = file
        self.method = method
        self.output = output
    
    def run_calculation(self):
        atoms = read( self.file )
        atoms.calc = TBLite( method = self.method )
        dynamics = BFGS( atoms, trajectory = self.output )
        dynamics.run( fmax = 0.05 )