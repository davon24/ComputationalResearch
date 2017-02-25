from numpy import *
import matplotlib.pyplot as plt


graphtype = "0005interp.txt"

data = loadtxt(graphtype, skiprows = 3) #file name with data points
[x,y] = [data[:,0], data[:,1]]
n = 4
i1 = 0
#initialize some variables
numsides = 4 #beginning with 8x8 grid
numtimes = 10 #each iteration changes numsidesx2 16x16 32x32 ect..

#initialize array
arrayx_y = zeros(numtimes)
boxsizes = zeros(numtimes)

#calculate the range of the data points for all three variables
rangex = max(x)-min(x)
rangey = max(y)-min(y)

#set the minimum value of the data points to a constant for all three variables
minx = min(x)
miny = min(y)

j=0
while j < numtimes:
    
    #initialize the box arrays for all three variables
    box_x_y = zeros([numsides,numsides])

    #equation to assign each data point to a region of the box array
    datasize = data.shape
    for i in range(datasize[0]):
        u = int(((x[i] - minx)/rangex)*(numsides))
        v = int(((y[i] - miny)/rangey)*(numsides))

    #correction for the maximum value in the data points   
        if u == numsides:
            u = u - 1
        if v == numsides:
            v = v - 1

        box_x_y[u,v] =1
        
    arrayx_y[j] = sum(box_x_y)
    boxsizes[j] = (numsides)
    numsides = numsides*2
    
    j = j+1
sumxy = 0
sumx = 0
sumy = 0
sumx2 = 0
for i in range(i1,n+i1):
    sumxy += log(boxsizes[i]) * log(arrayx_y[i])
    sumx += log(boxsizes[i])
    sumy += log(arrayx_y[i])
    sumx2 += (log(boxsizes[i]))**2

slope =( n*sumxy - sumx*sumy)/(n*sumx2 - (sumx)**2)
textstr = 'fractal=%.5f\n%s' %(slope, graphtype)

linex = array([boxsizes[i1],boxsizes[numtimes-1]])
liney = array([arrayx_y[i1], arrayx_y[i1]*(boxsizes[numtimes-1]/boxsizes[i1])**slope])



fig, ax = plt.subplots(1)
ax.loglog(boxsizes,arrayx_y,'bo')
ax.set_title ("Fractal box count")
ax.set_ylabel ("Number of occupied grid points")
ax.set_xlabel ("Number of grid points on a side")
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
ax.plot(linex,liney)
plt.show()
plt.figure(2)
plt.plot(x,y,'bo')
plt.show()

