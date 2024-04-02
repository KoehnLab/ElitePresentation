import py3Dmol
from ase.io.trajectory import Trajectory
from ase.io import read, write
import os

class TrajectoryVisualizer:
    
    view = py3Dmol.view(width=400, height=300)

    def __init__(self, file):
        self.file = file

    def create_trajectory(self):

        traj = Trajectory(self.file)

        #get current working directory and make a scratch 
        #directory
        path = os.getcwd()
        path = path + '/scratch'
        if not os.path.exists(path): os.makedirs(path)
        
        #output file name
        outFileName = 'trajectory.xyz'

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
            outFile = open(outFileName, 'a')
            outFile.write(fileStr)

    def show_trajectory(self):

        self.create_trajectory()

        xyz_data = open('trajectory.xyz',"r").read()

        self.view.addModelsAsFrames(xyz_data)
        self.view.animate({"loop": "forward"})
        # visualize with the stick option - can also consider spheres and more
        self.view.setStyle({"stick": {'colorscheme':'greyCarbon'}})
        self.view.show()