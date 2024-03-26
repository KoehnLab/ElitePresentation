from ase.io import read, write
from tblite.ase import TBLite
from ase.optimize import BFGS



def run_calculation(str: input, method = "GFN1-xTB", output = "trajectory.traj"):
    atoms = read( input )
    atoms.calc = TBLite( method = method )
    dynamics = BFGS( atoms, trajectory = output )
    dynamics.run( fmax = 0.05 )