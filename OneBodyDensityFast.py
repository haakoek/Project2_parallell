from numpy import *
import matplotlib.pyplot as plt

file_names = ["positionsIntN=10.txt","positionsNonN=10.txt"]

SIZE = 41
r_vec = linspace(0,4,SIZE)
step = r_vec[1] - r_vec[0]

for name in file_names:
	
	infile = open(name,'r')
	positions = []

	for line in infile:
		positions.append(float(line))

	Norm = len(positions)
	positions = array(positions)
	positions = sort(positions)

	max_r = max(positions)

	dr = 0.0
	prob = zeros(SIZE)

	counter = 0
	p = 0

	while p < len(positions):
	
		
		if(dr <= positions[p] < dr+step):
			prob[counter] += 1
			p += 1
		else:
			counter += 1
			dr += step
		    	 
	Norm = len(positions)
	prob = prob/(float(Norm))		

	plt.plot(r_vec,prob)

plt.legend(['Interacting','Non-interacting'])	
plt.show()

