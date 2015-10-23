__author__ = 'matthew'
import math
import Glob
import numpy as np


class Node(object):

    def __init__(self):

        # geometric variables
        self.coor = None
        self.isBound = False
        self.volu = None
        self.area = None

        # Variable for storing the boundary condition type
        # 1 = True and 0 = False
        self.boundType = {"Pressure": False, "Velocity": False}

        # solver variables
        self.velo = None
        self.press = None


class BoundFace(object):

    def __init__(self):
        self.bfNode = None
        self.bfCoef = None
        self.bfType = {"PressGrad": False, "PressValue": False , "VeloValue": False, "VeloGrad": False}


class Edge(object):

    def __init__(self):
        self.edgeCoef = None
        self.edgeLeng = None
        self.node0 = None
        self.node1 = None


class Mesh(object):

    def __init__(self):
        self.nodeTable = np.zeros(0, dtype=Node)
        self.edgeTable = np.zeros(0, dtype=Edge)
        self.bfTable = np.zeros(0, dtype=BoundFace)

        self.buildNodeTable()
        self.buildEdgeTable()
        self.buildBFTable()


    def buildNodeTable(self):

        self.nodeTable.resize(Glob.nx)
        xCoor = np.linspace(1.0, 1.0 + Glob.leng, Glob.nx)

        for inode in range(Glob.nx):
            self.nodeTable[inode] = Node()
            self.nodeTable[inode].coor = xCoor[inode]

        self.setBoundNodes()
        self.calcNodeVolume()
        self.calcNodeArea()

    def setBoundNodes(self):

        # First Node
        self.nodeTable[0].isBound = True
        self.nodeTable[0].boundType['Velocity'] = True

        # Last Node
        self.nodeTable[Glob.nx-1].isBound = True
        self.nodeTable[Glob.nx-1].boundType['Pressure'] = True

    def calcNodeVolume(self):

        for inode in range(Glob.nx):

            # Variables for the left and right x-values for node
            x1 = None
            x2 = None

            if inode == 0:
                x1 = self.nodeTable[inode].coor
                x2 = 0.5*(self.nodeTable[inode].coor + self.nodeTable[inode + 1].coor)
            elif inode == Glob.nx-1:
                x1 = 0.5*(self.nodeTable[inode].coor + self.nodeTable[inode -1].coor)
                x2 = self.nodeTable[inode].coor
            else:
                x1 = 0.5*(self.nodeTable[inode].coor + self.nodeTable[inode - 1].coor)
                x2 = 0.5*(self.nodeTable[inode].coor + self.nodeTable[inode + 1].coor)

            self.nodeTable[inode].volu = math.log(x2/x1)

    def calcNodeArea(self):

        for inode in range(Glob.nx):
            self.nodeTable[inode].area = 1.0/self.nodeTable[inode].coor

    def buildEdgeTable(self):

        self.edgeTable.resize(self.numEdges())

        for iedge in range(self.numEdges()):
            self.edgeTable[iedge] = Edge()
            self.edgeTable[iedge].node0 = iedge
            self.edgeTable[iedge].node1 = iedge + 1

        self.calcEdgeLeng()
        self.calcEdgeCoef()

    def calcEdgeLeng(self):

        for iedge in range(self.numEdges()):

            # Cache the nodes
            node0 = self.edgeTable[iedge].node0
            node1 = self.edgeTable[iedge].node1

            self.edgeTable[iedge].edgeLeng = self.nodeTable[node1].coor - self.nodeTable[node0].coor

    def numEdges(self):
        return Glob.nx - 1

    def calcEdgeCoef(self):

        for iedge in range(self.numEdges()):

            # Cache the nodes
            node0 = self.edgeTable[iedge].node0
            node1 = self.edgeTable[iedge].node1

            edgemidpoint = self.nodeTable[node0].coor + 0.5*self.edgeTable[iedge].edgeLeng
            edgeArea = 1.0/edgemidpoint

            self.edgeTable[iedge].edgeCoef = edgeArea

    def buildBFTable(self):
        self.bfTable.resize(2)

        # Inlet BF
        self.bfTable[0] = BoundFace()
        self.bfTable[0].bfNode = 0
        self.bfTable[0].bfCoef = -self.nodeTable[0].area
        self.bfTable[0].bfType['VeloValue'] = True

        # Outlet BF
        self.bfTable[1] = BoundFace()
        self.bfTable[1].bfNode = Glob.nx-1
        self.bfTable[1].bfCoef = self.nodeTable[Glob.nx-1].area
        self.bfTable[1].bfType['PressValue'] = True

