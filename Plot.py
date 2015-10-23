__author__ = 'matthew'
import numpy as np
from Glob import nx as nx
import matplotlib.pyplot as plt

class Plotter(object):

    def __init__(self, mesh):

        # Variable holders
        self.velo =  np.zeros(nx, dtype=float)
        self.press = np.zeros(nx, dtype=float)
        self.x = np.zeros(nx, dtype=float)
        self.energy = np.zeros(nx, dtype=float)

        self.getVariable(mesh)

    def getVariable(self, mesh):

        for inode in range(nx):
            self.x[inode] = mesh.nodeTable[inode].coor
            self.velo[inode] = mesh.nodeTable[inode].velo
            self.press[inode] = mesh.nodeTable[inode].press

            self.energy[inode] = self.press[inode] + 0.5*self.velo[inode]*self.velo[inode]

    def plot(self):
        plt.plot(self.x, self.press, label="Pressure")
        plt.plot(self.x, self.velo, label="Velo")
        plt.plot(self.x, self.energy, label="Energy")

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
        plt.show()

