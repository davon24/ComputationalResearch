from numpy import *
from percolationfunctions2 import *
import matplotlib.pyplot as plt
import numpy.random as random
import time as time





def MSTfunction(points, L,avgpointsperbox=3):
#This function takes the following inputs:
#points: A list of coordinate pairs (x,y) as floating point arrays
#L: For wrapping purposes, the size of the window that these points sit in.
#   We are working in a toroidal window that wraps around
#avgpointsperbox is a target number of points to put in a box.  We will have AT MOST
#   that many points in a box; we may have fewer.
#This function is optimized for large numbers of points, so that we only compute distances
#to points in their immediate neighborhood


    #Start by defining the list that will hold the mst.  It will become a list of links,
    #The link in each row being the link from the corresponding point to the mst
    N = len(points)
    mst = [[] for i in range(N)]
    mst[0] = link(0, 0, array([0.0,0.0]), array([0.0,0.0]))
    #How big should the box be?
    w = int(sqrt(N/avgpointsperbox))
    boxwidth = L/w
    
    #Create a box full of boxes
    boxes = [[ [] for i in range(w)] for j in range(w)]
    #sort points into boxes
    for k in range(1,N):
        x = int(points[k][0]/boxwidth)
        y = int(points[k][1]/boxwidth)
        boxes[x][y]+=[k]

    #This is the list of distances from points not on MST to points on the MST
    distances = []
    points_on = []
    points_off = []
    
    for i in range(1,N): # finds the distances and links between point 0 and the rest
        points_on = points_on + [0]
        points_off = points_off + [i]
        distances.append([wrapping_squared_dist(points[0][0],points[0][1],points[i][0],points[i][1],L)])
 
    #We have to add another N-1 points to the MST
    min_finding_time=0
    rest_of_time = 0
    for t in range(N-1):
        t1=time.time()
        mindist = min(distances)
        t2=time.time()
        min_finding_time=min_finding_time+t2-t1
        #Where is the shortest distance?
        row = distances.index(mindist)
        #Figure out the point we're adding to the MST
        new_point = points_off[row]
        #Figure out the existing point that it's connecting to
        old_point = points_on[row]
        xold = points[old_point][0]
        xnew = points[new_point][0]
        yold = points[old_point][1]
        ynew = points[new_point][1]
	
	#Now we compute the bond vector, taking into account the possibility of wrapping
        if xnew - xold >= L/2:
            delta_x = xnew - xold - L
        elif xnew - xold < -L/2:
            delta_x = xnew - xold + L
        else:
            delta_x = xnew - xold
 
        if ynew - yold >= L/2:
            delta_y = ynew - yold - L
        elif ynew - yold < -L/2:
            delta_y = ynew - yold + L
        else:
            delta_y = ynew - yold
 
        mst[new_point] = link(new_point, old_point, array([delta_x, delta_y]), mst[old_point].rootvector + array([delta_x, delta_y])) 
    	#Now remove the point that we just added
        del distances[row] 
        del points_on[row]
        del points_off[row]
        #Remove that point from its box
        i = int(xnew/boxwidth)
        j = int(ynew/boxwidth)
        box=boxes[i][j]
        k = box.index(new_point)
        del box[k]
        boxes[i][j]=box
        
 
        #Compute distances from points not yet on MST to new point added to MST
        xvals = [(i-1)%w, i, (i+1)%w]
        yvals = [(j-1)%w, j, (j+1)%w]
        for x in xvals:
            for y in yvals:
                box = boxes[x][y]
                for point_off in box:
                    d = wrapping_squared_dist(points[new_point][0],points[new_point][1],points[point_off][0],points[point_off][1],L)
                    #Next we need to find the corresponding row in "distances", "points_on", and "points_off"
                    k = points_off.index(point_off)
                    if d<distances[k]:
                        distances[k] = d
                        points_on[k] = new_point
        t3=time.time()
        rest_of_time = rest_of_time + t3-t2
        
    #print('We found the minimum distance in the list %d times at an average of %.2E seconds per iteration.' % (N, min_finding_time/N))
    #print('The time for the rest of the stuff in the loop averaged %.2E seconds per iteration.' % (rest_of_time/N))
    return mst     

# starts creating the graph
numpoints = 1000 # Number of points on graph
numsets = 100 # Number of MSTs created
N = 100 # Number of points on each MST
L = 10.0 #This is the size of our toroidal window for wrapping purposes
random.seed(42)
critarray = zeros(numsets) # array where the critical points will be stored
perarray = zeros(numpoints) # array that holds the percentages

for i in range(numsets):
    points = [L*random.rand(2) for m in range(N)]  
    newMST = MSTfunction(points,L)
    

    #Find the extremes
    [pointA,pointB] = graph_extrema(newMST)
    #print('The extreme points are %d and %d.' %(pointA, pointB))
    #Now we find the path between extremes
    extrema_links = path_between_extremes(newMST,pointA,pointB)
    #Change 5,6 to 1,7 to get the left-right path  instead of top-bottom
    bondlengths=[A.bondlengthsquared for A in extrema_links]
    criticalbond = sqrt(max(bondlengths))
    critarray[i] = criticalbond

steps = amax(critarray)/numpoints
grapharray = ([i*steps for i in range(numpoints)])

for i in range(numpoints):
    for j in range(numsets):
        if critarray[j] <= grapharray[i]:
            perarray[i] += 1.0

perarray /= numsets

plt.plot(grapharray, perarray, 'bo')
plt.xlabel("Critical bond length")
plt.ylabel("Percent")
plt.show()
print critarray
print perarray