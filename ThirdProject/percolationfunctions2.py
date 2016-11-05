from numpy import *
import matplotlib.pyplot as plt

#Function for finding distance in a square 2D wrapping domain

class link: #Used to represent links to other sites in the minimal spanning tree
    def __init__(self,start,parent=[],bondvector=array([0.0,0.0]),rootvector=array([0.0,0.0])):
    #Start is the index of the site joining the mst
    #Default parent is empty, for the case of a root site
    #Default bondvector is zero, for cases where we don't need it.
        self.start = start
        self.parent = parent
        self.bondvector = bondvector
        self.rootvector = rootvector
        self.bondlengthsquared = dot(bondvector,bondvector)
        self.rootlengthsquared = dot(rootvector,rootvector)

def wrapping_squared_dist(x1,y1,x2,y2,L=0):
    #We assign L the default value 1, so that if it isn't specified we assume
    #We are working in the unit square
    dysquare = min([(y1-y2)**2, (min([y1,y2])+L-max([y1,y2]))**2])
    dxsquare = min([(x1-x2)**2, (min([x1,x2])+L-max([x1,x2]))**2])
    return dxsquare + dysquare
    
def graph_MST(MST, points, L):
    n = len(points)
    graphx = [[points[i][0]] for i in range(n)]
    graphy = [[points[i][1]] for i in range(n)]
    plt.figure()
    ax = plt.axes()
    plt.axis([0,L,0,L])
    plt.plot(graphx, graphy, 'ro')

    for link in MST:
        [x1,y1] = points[link.parent]
        [x2,y2] = points[link.start]

        if (x2 - x1 >= L/2.0):
            delta_x = x2 - x1 - L
        elif x2 - x1 < -L/2.0:
            delta_x = x2 - x1 + L
        else:
            delta_x = x2 - x1

        if y2 - y1 >= L/2.0:
            delta_y = y2 - y1 - L
        elif y2 - y1 < -L/2.0:
            delta_y = y2 - y1 + L
        else:
            delta_y = y2 - y1

        
        plt.arrow(x1, y1, delta_x, delta_y, head_width=0.2, head_length=0.3, fc='k', ec='k')
        plt.arrow(x2, y2, -delta_x, -delta_y)

    for i in range(n):
        plt.annotate('%d' % (i), xy=points[i], xycoords='data',
                xytext=(-3, 5), textcoords='offset points'
                )
    return plt.show()
 
def graph_extrema(mst):
#December 4, 2015, written by Alex Small
#Takes a list of links and returns the two points that are the farthest-separated along either the x or y direction
    xvalues = [L.rootvector[0] for L in mst]
    yvalues = [L.rootvector[1] for L in mst]
    
    dx = max(xvalues)-min(xvalues)
    dy = max(yvalues)-min(yvalues)
    
    if dx>dy:
        pointA = xvalues.index(max(xvalues))
        pointB = xvalues.index(min(xvalues))
    else:
        pointA = yvalues.index(max(yvalues))
        pointB = yvalues.index(min(yvalues))
    
    return [pointA, pointB]

def path_between_extremes(mst, point1, point2, root=0):
    #November 24, 2015, written by Alex Small
    #This function finds a path between the two extremes of the MST
    #Inputs:
    #mst = list of link objects, ordered such that the link from point i to the
    #       mst is in row i of the list
    #point 1 = index of one of the extreme points
    #point 2 = index of the other extreme point
    #root = optional index of root site (assumed to be 1 unless otherwise specified)


    #find the links from parent 1 to root
    from1toroot = [mst[point1]] 
    parent = mst[point1].parent
    while parent != root:
        from1toroot = [mst[parent]]+from1toroot #Add the next link to the mst
        parent = mst[parent].parent #Find the link between the next site and the mst
    
    #Now do the same for point2
    from2toroot = [mst[point2]] #Start with an empty list of links
    parent = mst[point2].parent
    while parent != 0:
        from2toroot = [mst[parent]]+from2toroot #Add the next link to the mst
        parent = mst[parent].parent #Find the link between the next site and the mst

    #Now go through the lists of links and remove identical elements
    while from1toroot[0]==from2toroot[0]:
        # print('First bond in from1toroot joins %d and %d.\n' % (from1toroot[0].start, from1toroot[0].parent))
        # print('First bond in from2toroot joins %d and %d.\n' % (from2toroot[0].start, from2toroot[0].parent))

        if len(from1toroot)==1 or len(from2toroot)==1:
            break

        del from1toroot[0]
        del from2toroot[0]
        # print('We just did a deletion.\n')

    
    return from1toroot+from2toroot
        
        
