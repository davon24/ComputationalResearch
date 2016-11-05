from numpy import *
from scipy import *
from matplotlib.pyplot import *
# from pylab import * # loads matplotlib

N = 20
A_1,A_2 = meshgrid( arange(0,N,1),arange(0,N,1) )
rho = zeros((N,N)) # potential density
U = zeros((N,N))
E_x = zeros((N,N))
E_y = zeros((N,N))

rho[7,10]  =  10


numcharge = 2
maxdiff = 0.000000000001
testdiff = 1
diff = 0
U = rho*1.0

def interpolate(E, x, y):
	delta = x-floor(x)
	beta  = y-floor(y)
	I1 = (delta-1.0)*E[int(x),int(y)+1] + (delta)*E[int(x)+1,int(y)+1]
	I2 = (delta-1.0)*E[int(x),int(y)] + (delta)*E[int(x)+1,int(y)]
	I3 = (beta-1.0)*I1 + (beta)*I2
	return I3

while testdiff > maxdiff:
	testdiff = 0.0
	for i in range(1,N-1):
		for j in range(1,N-1):
			U_0 = (U[i-1,j]+U[i+1,j]+U[i,j-1]+U[i,j+1])/4+rho[i,j]/4
			diff  = abs(U_0-U[i,j])
			if diff > testdiff:
				testdiff = diff
			U[i,j] = U_0
for i in range(1,N-1):
	for j in range(1,N-1):
		E_y[i,j] = -(-(U[i+1,j]-U[i-1,j])/2.0)    # the first index is the y-axis
		E_x[i,j] = -(-(U[i,j+1]-U[i,j-1])/2.0)   # second index is the x-axis

# interpolate the electric field at point (x,y)
# x = 9.5
# y = 9.5
# E_t = zeros((2))
# E_t[0] = interpolate(E_x , x, y)
# E_t[1] = interpolate(E_y , x, y)

# plotting function
plt=imshow(U,interpolation = "nearest")
xlim([19,0])
ylim([19,0])
title('Electric Potential')
# add a colorbar for pseudoprobability
cbar=colorbar(plt)#create colorbar
cbar.ax.set_ylabel('Electric Potential')
xlabel('x')
ylabel('y')
# show plot
show()
print U
# Q = quiver( A_1, A_2, E_x, E_y,
#             pivot='mid', color='r', scale=1/0.50 )
# plt = quiverkey(Q, 0.5, 0.03, 0, r'', fontproperties={'weight': 'bold'})
# plot( A_1, A_2, 'k.')
# xlim([19,0])
# ylim([19,0])
# title("Electric Field")
# show()