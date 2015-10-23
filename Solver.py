__author__ = 'matthew'
from Mesh import Mesh
import Glob
import numpy as np


class Solver:

    def __init__(self, mesh):

        # Data members
        self.AMatr = np.zeros((Glob.nx, Glob.nx),dtype=float)
        self.bVect = np.zeros(Glob.nx, dtype=float)

        self.veloGrad = np.zeros(Glob.nx, dtype=float)
        self.pressGrad = np.zeros(Glob.nx, dtype=float)

        self.press = np.zeros(Glob.nx, dtype=float)
        self.deltaWStar = np.zeros(Glob.nx, dtype=float)

    def setInitCond(self, mesh):
        for inode in range(Glob.nx):

            # Setting initial BC's
            if mesh.nodeTable[inode].isBound:
                if mesh.nodeTable[inode].boundType["Pressure"]:
                    mesh.nodeTable[inode].press = Glob.P_out
                    mesh.nodeTable[inode].velo = 0.0
                elif mesh.nodeTable[inode].boundType["Velocity"]:
                    mesh.nodeTable[inode].velo = Glob.U_in
                    mesh.nodeTable[inode].press = 0.0
            else:
                mesh.nodeTable[inode].press = 0.0
                mesh.nodeTable[inode].velo = 0.0

            # mesh.nodeTable[inode].press = 2.0 - 0.5*mesh.nodeTable[inode].coor*mesh.nodeTable[inode].coor
            # mesh.nodeTable[inode].velo = mesh.nodeTable[inode].coor

    def stepZero(self, mesh):

        # Zero the vector
        self.deltaWStar[:] = 0.0


        for iedge in range(mesh.numEdges()):

            # Cache the nodes
            node0 = mesh.edgeTable[iedge].node0
            node1 = mesh.edgeTable[iedge].node1

            faceFlux = 0.0

            veloUpwind = mesh.nodeTable[node0].velo
            veloFace = 0.5*(mesh.nodeTable[node0].velo + mesh.nodeTable[node1].velo)
            faceFlux = -Glob.deltaT*veloUpwind*veloFace*mesh.edgeTable[iedge].edgeCoef

            self.deltaWStar[node0] += faceFlux
            self.deltaWStar[node1] -= faceFlux

        for bf in range(mesh.bfTable.size):

            #bfnode
            bfnode = mesh.bfTable[bf].bfNode

            veloUpwind = mesh.nodeTable[bfnode].velo
            veloFace = mesh.nodeTable[bfnode].velo
            faceFlux = -Glob.deltaT*veloUpwind*veloFace*mesh.bfTable[bf].bfCoef

            self.deltaWStar[bfnode] += faceFlux

        for inode in range(Glob.nx):
            self.deltaWStar[inode] /= mesh.nodeTable[inode].volu

    def stepOne(self, mesh):

        #Zero the a and b matricies
        self.AMatr[:, :] = 0.0
        self.bVect[:] = 0.0

        for iedge in range(mesh.numEdges()):

            # Cache the nodes
            node0 = mesh.edgeTable[iedge].node0
            node1 = mesh.edgeTable[iedge].node1

            #variable for faceflux
            faceflux = 0.0

            # div(A*U_x)
            uface = 0.5*(mesh.nodeTable[node0].velo + mesh.nodeTable[node1].velo)
            faceflux += (uface*mesh.edgeTable[iedge].edgeCoef)

            # div(A*U_x*d/dx(u_x))
            gradU = (mesh.nodeTable[node1].velo - mesh.nodeTable[node0].velo)/mesh.edgeTable[iedge].edgeLeng
            faceflux -= Glob.deltaT*uface*gradU*mesh.edgeTable[iedge].edgeCoef

            #deltaWStarFace
            # deltaWFace = 0.5*(self.deltaWStar[node0] + self.deltaWStar[node1])
            # faceflux += deltaWFace*mesh.edgeTable[iedge].edgeCoef

            # add to nodes
            self.bVect[node0] += faceflux
            self.bVect[node1] -= faceflux

            # edge Press coef
            pressCoef = mesh.edgeTable[iedge].edgeCoef/mesh.edgeTable[iedge].edgeLeng

            # add to AMatr
            self.AMatr[node0][node0] -= pressCoef
            self.AMatr[node0][node1] += pressCoef

            self.AMatr[node1][node1] -= pressCoef
            self.AMatr[node1][node0] += pressCoef

        # Multiple pressure matrix by deltaT
        self.AMatr *= Glob.deltaT

        # Account for boundary condition
        for bf in range(mesh.bfTable.size):

            # cache the BF node
            bfnode = mesh.bfTable[bf].bfNode

            if mesh.bfTable[bf].bfType["PressValue"]:
                self.AMatr[bfnode,:] = 0.0
                self.AMatr[bfnode,bfnode] = 1.0
                self.bVect[bfnode] = Glob.P_out

            if mesh.bfTable[bf].bfType["VeloValue"]:
                faceVelo = mesh.nodeTable[bfnode].velo
                self.bVect[bfnode] += faceVelo*mesh.bfTable[bf].bfCoef

        self.press = np.linalg.solve(self.AMatr,self.bVect)

        for inode in range(Glob.nx):
            mesh.nodeTable[inode].press = self.press[inode]

    def calcPressGrad(self, mesh):

        #zero pressGrad vectors first
        self.pressGrad[:] = 0.0

        for iedge in range(mesh.numEdges()):

            # cache nodes
            node0 = mesh.edgeTable[iedge].node0
            node1 = mesh.edgeTable[iedge].node1

            pressFace = 0.5*(mesh.nodeTable[node0].press + mesh.nodeTable[node1].press)

            self.pressGrad[node0] += pressFace#*mesh.edgeTable[iedge].edgeCoef
            self.pressGrad[node1] -= pressFace#*mesh.edgeTable[iedge].edgeCoef

        for bf in range(mesh.bfTable.size):

            #bf node
            bfnode = mesh.bfTable[bf].bfNode

            pressFace = mesh.nodeTable[bfnode].press

            if bf is 0:
                self.pressGrad[bfnode] -= pressFace
            elif bf is 1:
                self.pressGrad[bfnode] += pressFace

        for inode in range(Glob.nx):
            self.pressGrad[inode] /= mesh.nodeTable[inode].volu

    def calcVeloGrad(self, mesh):

        # First zero vector
        self.veloGrad[:] = 0.0

        for iedge in range(mesh.numEdges()):

            # cache nodes
            node0 = mesh.edgeTable[iedge].node0
            node1 = mesh.edgeTable[iedge].node1

            veloFace = 0.5*(mesh.nodeTable[node0].velo + mesh.nodeTable[node1].velo)

            self.veloGrad[node0] += veloFace
            self.veloGrad[node1] -= veloFace

        for bf in range(mesh.bfTable.size):

            #cache bfNode
            bfnode = mesh.bfTable[bf].bfNode

            veloFace = mesh.nodeTable[bfnode].velo

            if bf is 0:
                self.veloGrad[bfnode] -= veloFace
            elif bf is 1:
                self.veloGrad[bfnode] += veloFace

        for inode in range(Glob.nx):
            self.veloGrad[inode] /= mesh.nodeTable[inode].volu

    def stepTwo(self, mesh):

        for inode in range(Glob.nx):
            if mesh.nodeTable[inode].boundType["Velocity"] is False:
                mesh.nodeTable[inode].velo = mesh.nodeTable[inode].velo - Glob.deltaT*self.pressGrad[inode] \
                                            -mesh.nodeTable[inode].velo*Glob.deltaT*self.veloGrad[inode]

    # mass cosveration d/dx(A*U_x) = 0
    def checkMassConsv(self, mesh):

        massImbal = np.zeros(Glob.nx, dtype=float)

        for iedge in range(mesh.numEdges()):

            # cache nodes
            node0 = mesh.edgeTable[iedge].node0
            node1 = mesh.edgeTable[iedge].node1

            veloFace = 0.5*(mesh.nodeTable[node0].velo + mesh.nodeTable[node1].velo)

            massImbal[node0] += veloFace*mesh.edgeTable[iedge].edgeCoef
            massImbal[node1] -= veloFace*mesh.edgeTable[iedge].edgeCoef

        for bf in range(mesh.bfTable.size):

            #bfNode
            bfnode = mesh.bfTable[bf].bfNode

            bfVelo = mesh.nodeTable[bfnode].velo

            massImbal[bfnode] += bfVelo*mesh.bfTable[bf].bfCoef

        for inode in range(Glob.nx):
            massImbal[inode] /= mesh.nodeTable[inode].volu

        return np.linalg.norm(massImbal)

    def solverLoop(self, mesh):

        self.setInitCond(mesh)

        for iter in range(Glob.iter):
            # self.stepZero(mesh)
            self.stepOne(mesh)
            self.calcPressGrad(mesh)
            self.calcVeloGrad(mesh)
            self.stepTwo(mesh)
            print ("Iter:\t" + str(iter) + "\t" "Mass Imbal: " + str(self.checkMassConsv(mesh)))


