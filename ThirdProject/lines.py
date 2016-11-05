from numpy import *
from percolationfunctions2 import *
import matplotlib.pyplot as plt
import random
import time as time

random.seed(3)
L = 4 #length and width of window
N = 5 #number of lines
P = 6 #number of points

starting = zeros([N,3])
ending = zeros([N,3])
points = zeros([P,3])

for i in range(N):	#creates the starting points of the lines
	for j in range(2):
		starting[i,j] = L*random.random() 
	starting[i,2] = pi*random.triangular(-1,1)

for i in range(N):	#ending points of the lines
	theta = starting[i,2]
	ending[i,0] = starting[i,0] + L*cos(theta)
	ending[i,1] = starting[i,1] + L*sin(theta)
	ending[i,2] = starting[i,2]

for i in range(P):	#creates the points that will go on the lines
	s = random.random()
	index = random.randint(0,N-1)
	xp = s*starting[index,0] + (1-s)*ending[index,0]
	yp = s*starting[index,1] + (1-s)*ending[index,1]
	points[i,0] = xp
	points[i,1] = yp  
	points[i,2] = index


print starting
print ending
print sqrt((starting[0,0]-ending[0,0])**2 + (starting[0,1]-ending[0,1])**2) 
print points

for i in range(N):
	x = array([starting[i,0],ending[i,0]])
	y = array([starting[i,1],ending[i,1]])
	plt.plot(x,y,'b')

plt.plot(points[:,0],points[:,1],'ro')

plt.show()