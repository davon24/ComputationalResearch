import math
import numpy as np
from scipy.special import jn
from matplotlib.pyplot import *
import numpy.random.mtrand as mt
import numpy.random as nr

nr.seed(6)
poisson = mt.poisson
x0 = 3   # x coordinate of the origin
y0 = 2   # y coordinate of the origin
wavelength = 6    # wavelength in pixels
photons = 500

formulameans = np.zeros([7,7]) # creates an array of zeros

def bessel(x0,y0,wavelength):             # defines our formula as bessel
    k = (math.pi*2)/wavelength
    r = math.sqrt(((x2-x0)**2) + ((y2-y0)**2))
    p = ((2*jn(1,(k*r))/(k*r))**2)
    return p


for b in range(7):              #creates a 7x7 grid composed of 4,900 subpixels 
    for a in range(7):
        rsum = 0
        for y in range(10):
            y2 = (((0.1*y)-0.45) + (1*b))
            for x in range(10):
                x2 = (((0.1*x)-0.45) + (1*a))
                rsum = (rsum + bessel(x0,y0,wavelength))
        average = ((rsum)/100)
        formulameans[b,a] = average  

Isum = formulameans.sum()           # adding all the contents of the array called formulameans
picture = formulameans/Isum         # dividing everything in formulameans by the sum and creating an array named picture
normalized = picture*photons            # multiplying picture by the desired number of photon counts 
normalized += 0                         # adding a number to the photon count of each pixel
graph = poisson(normalized)             # replaces the array of photon counts with an array of poisson random numbers


plt=imshow(normalized,interpolation = "nearest")  # creates the color-coded graph
xlim([0,6])
ylim([0,6])
title('')
cbar=colorbar(plt)
cbar.ax.set_ylabel('')
xlabel('')
ylabel('')
show()

plt=imshow(graph,interpolation = "nearest")
xlim([0,6])
ylim([0,6])
title('')
cbar=colorbar(plt)
cbar.ax.set_ylabel('')
xlabel('')
ylabel('')
show()



xsum = 0.0
x2sum = 0.0
graphcount = 1000     # number of graphs that are being created

for i in range(graphcount):
    cen = 0.0
    image = poisson(normalized)
    imagesum = image.sum()
    for y in range(7):
        for x in range(7):
            cen = cen + image[y,x] * (x+1)
    xcen = cen/imagesum
    xsum = xsum + xcen
    x2sum = x2sum + xcen**2
    
xavg = xsum/graphcount
x2avg = x2sum/graphcount

print xavg
print x2avg

standarddeviation = math.sqrt(x2avg - (xavg)**2)
print standarddeviation