from math import*
from numpy import*
from scipy.special import jn
from matplotlib.pyplot import*
import numpy.random.mtrand as mt
import numpy.random as nr
nr.seed(3)
poisson = mt.poisson
x0 = 4  # x coordinate of the origin
y0 = 4  # y coordinate of the origin
wavelength = 6.0    # wavelength in pixels
photons = 2000
ky = (2.0*pi)/wavelength
kx = 1.5*pi/wavelength
formulameans = zeros([9,9]) # creates an array of zeros
b = 2.0

def airyPSF(x, y,wavelength):             # defines the diffraction-limited PSF
#Inputs are relative coordinates, coordinates relative to molecular position
    kr = sqrt(((kx**2)*(x**2)) + ((ky**2)*(y**2)))
    p = ((2*jn(1,(kr))/(kr))**2)
    return p

numsubpix = 10
for c in range(9):              #creates a 7x7 grid composed of 4,900 subpixels 
    for a in range(9):
        rsum = 0.0
        for y in range(numsubpix):
            y2 = (((0.1*y)-0.45) + (1*c))
            for x in range(numsubpix):
                x2 = (((0.1*x)-0.45) + (1*a))
                rsum = (rsum + airyPSF(x2-x0,y2-y0,wavelength))
        avgcount = ((rsum)/numsubpix**2)
        formulameans[c,a] = avgcount

Isum = formulameans.sum()           # adding all the contents of the array called formulameans
picture = formulameans/Isum         # dividing everything in formulameans by the sum and creating an array named picture
normalized = picture*photons
normalized += b           # multiplying picture by the desired number of photon counts 
                         # adding a number to the photon count of each pixel
S = poisson(normalized)             # replaces the array of photon counts with an array of poisson random numbers
I0 = S.max()

plt=imshow(formulameans,interpolation = "nearest")
xlim([0,8])
ylim([0,8])
title('')
cbar=colorbar(plt)
cbar.ax.set_ylabel('')
xlabel('')
ylabel('')
show()


def psf(x0,y0,I0,b,kx,ky,x,y):
    ans = b + (I0 * (e**(-0.25*(((kx**2)*((x-x0)**2)) + ((ky**2)*((y-y0)**2))))))
    return ans

  
diff1 = 1
diff2 = 1

def dx(x,y): 
    ans = I0*0.5*(kx**2)*(x-x0)*(e**(-0.25*(((kx**2)*((x-x0)**2)) + ((ky**2)*((y-y0)**2)))))
    return ans
def dy(x,y):
    ans = I0*0.5*(ky**2)*(y-y0)*(e**(-0.25*(((kx**2)*((x-x0)**2)) + ((ky**2)*((y-y0)**2)))))
    return ans
def dx2(x,y):
    ans = I0 * ((0.25*(kx**4)*((x-x0)**2)) - (0.5*(kx**2))) * (e**(-0.25*(((kx**2)*((x-x0)**2)) + ((ky**2)*((y-y0)**2)))))
    return ans                                                                                                   # Defining the derivatives of the loglikelihood function
def dy2(x,y):
    ans = I0 * ((0.25*(ky**4)*((y-y0)**2)) - (0.5*(ky**2))) * (e**(-0.25*(((kx**2)*((x-x0)**2)) + ((ky**2)*((y-y0)**2)))))
    return ans    
def dxy(x,y):
    ans = I0 * 0.25 * (kx**2) * (ky**2) * (x-x0) * (y-y0) * (e**(-0.25*(((kx**2)*((x-x0)**2)) + ((ky**2)*((y-y0)**2)))))
    return ans




# Finding the standard deviation    
xsum = 0.0
x2sum = 0.0
ysum = 0.0
y2sum = 0.0
graphcount = 1000     # number of graphs that are being created
old = np.matrix([[3],[3]])      # our initial guess for the pixel


for i in range(graphcount):
    
    S = poisson(normalized)        
    I0 = S.max()
    DX = 0.0
    DY = 0.0
    DX2 = 0.0
    DXY = 0.0
    DY2 = 0.0
    diff1 = 1
    diff2 = 1
    while diff1 > 0.05 or diff2 > 0.05:             
        x0 = old[0,0]
        y0 = old[1,0]
        for y in range(9):
           for x in range(9):
                DX += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y))-1)*dx(x,y)) 
                DY += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y))-1)*dy(x,y))
                DX2 += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y))-1)*dx2(x,y)) - (((dx(x,y))**2)*S[y,x]/((psf(x0,y0,I0,b,kx,ky,x,y))**2))
                DXY += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y))-1)*dxy(x,y)) - (dx(x,y)*dy(x,y)*S[y,x]/(psf(x0,y0,I0,b,kx,ky,x,y))**2)
                DY2 += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y))-1)*dy2(x,y)) - (((dy(x,y))**2)*S[y,x]/((psf(x0,y0,I0,b,kx,ky,x,y))**2))

        J = np.matrix([[DX2,DXY],[DXY,DY2]])
        J2 = J.getI()
        F = np.matrix([[DX],[DY]])
        F = -F
        D = J2*F
        new = old + D
        diff1 = abs(old[0,0] - new[0,0])
        diff2 = abs(old[1,0] - new[1,0])
        old = new
    print old, i
    xsum = xsum + old[0,0]
    x2sum = x2sum + old[0,0]**2
    ysum = ysum + old[1,0]
    y2sum = y2sum + old[1,0]**2
    
xavg = xsum/graphcount
x2avg = x2sum/graphcount
yavg = ysum/graphcount
y2avg = y2sum/graphcount
print xavg
#print x2avg
print yavg
#print y2avg

standarddeviationx = sqrt(x2avg - (xavg)**2)
standarddeviationy = sqrt(y2avg - (yavg)**2)

print standarddeviationx , standarddeviationy