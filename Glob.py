__author__ = 'matthew'
# Solver Global Variables
nx = 3
leng = 1.0
iter = 1

# Varibale data
U_in = 1.0
P_out = 0.0

#Calc deltaT
# C = U*dt/dx
CFL = 0.1
deltaT = CFL*(leng/(nx-1))/U_in


print "Starting Simulation"
print "DeltaT: ", deltaT, "\n"