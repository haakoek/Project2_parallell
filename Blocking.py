from numpy import *
import matplotlib.pyplot as plt

infile = open("energies.txt",'r')
local_energies = []

for line in infile:
	local_energies.append(float(line))


local_energies = array(local_energies)
N = len(local_energies)

M = []

for i in range(1,3000):
	
	if N%i == 0:
		M.append(N/i)
	
n_vec = zeros(len(M))
var_blocking = zeros(len(M))
k = 0

for m in M:
	
	n = N/m
	E_new = zeros(m)
	
	for j in range(0,m):
		E_new[j] = sum(local_energies[j*n:(j+1)*n])/n
		
	expectation_Enew = sum(E_new)/float(m)
	expectation_Enew2 = sum(E_new**2)/float(m)
	
	var_blocking[k] = (expectation_Enew2 - expectation_Enew**2)/m
	n_vec[k] = n
	k += 1
	
	print m
	
print var_blocking[0]*M[0]


plt.plot(n_vec,var_blocking,'-o')
plt.xlabel('N/m')
plt.ylabel('Var')
plt.show()

