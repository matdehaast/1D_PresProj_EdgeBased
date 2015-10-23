__author__ = 'matthew'
from Mesh import Mesh
from Solver import Solver
from Plot import Plotter
import Glob

#Create mesh
mesh = Mesh()

#create solver and run loop
solver = Solver(mesh=mesh)
solver.solverLoop(mesh)

#Create plotter and plot
plotter = Plotter(mesh)
plotter.plot()
