from math import*
from numpy import*
from scipy.special import jn
from matplotlib.pyplot import*
import numpy.random.mtrand as mt
import numpy.random as nr
nr.seed(3)
poisson = mt.poisson
x0 = 10  # x coordinate of the origin
y0 = 9  # y coordinate of the origin
wavelength = 4.0    # wavelength in pixels
photons = 1000000
kx = 1.0*pi/wavelength
ky = 1.0*pi/wavelength
imagesize = 19
formulameans = zeros([imagesize,imagesize]) # creates an array of zeros
b = 5
theta = 0
def PSF(x, y,wavelength):             # defines the diffraction-limited PSF
#Inputs are relative coordinates, coordinates relative to molecular position
	kr = sqrt((kx**2 * cos(theta)**2 + ky**2 * sin(theta)**2)*(x)**2 + (ky**2 * cos(theta)**2 + kx**2 * sin(theta)**2)*(y)**2 + (kx**2 - ky**2)*(x)*(y)*sin(2*theta))
	return e**(-0.25*(kr**2))

numsubpix = 10
for c in range(imagesize):              #creates a grid composed of subpixels 
    for a in range(imagesize):
        rsum = 0.0
        for y in range(numsubpix):
            y2 = (((0.1*y)-0.45) + (1*c))
            for x in range(numsubpix):
                x2 = (((0.1*x)-0.45) + (1*a))
                rsum = (rsum + PSF(x2-x0, y2-y0,wavelength) *((x2-x0)*2*pi/wavelength)**2)
        avgcount =  ((rsum)/numsubpix**2)+b
        formulameans[c,a] = avgcount

Isum = formulameans.sum()           # adding all the contents of the array called formulameans
picture = formulameans/Isum         # dividing everything in formulameans by the sum and creating an array named picture
normalized = picture*photons
normalized += b           # multiplying picture by the desired number of photon counts 
S = poisson(normalized)                   # adding a number to the photon count of each pixel
             # replaces the array of photon counts with an array of poisson random numbers



plt=imshow(S,interpolation = "nearest")
xlim([0,imagesize - 1])
ylim([0,imagesize - 1])
title('')
cbar=colorbar(plt)
cbar.ax.set_ylabel('')
xlabel('')
ylabel('')
show()

plt=imshow(normalized,interpolation = "nearest")
xlim([0,imagesize - 1])
ylim([0,imagesize - 1])
title('')
cbar=colorbar(plt)
cbar.ax.set_ylabel('')
xlabel('')
ylabel('')
show()

def beta(x0,y0,I0,b,kx,ky,x,y,theta):
	ans = e**(-0.25*((kx**2 * cos(theta)**2 + ky**2 * sin(theta)**2)*(x-x0)**2 + (ky**2 * cos(theta)**2 + kx**2 * sin(theta)**2)*(y-y0)**2 + (kx**2 - ky**2)*(x-x0)*(y-y0)*sin(2.0*theta)))
	return ans

def psf(x0,y0,I0,b,kx,ky,x,y,theta):
	ans = b + ((x-x0)*(2.0*pi/wavelength))**2 * I0*beta(x0,y0,I0,b,kx,ky,x,y,theta)
	return ans

  
diff1 = 1
diff2 = 1

def dx(x,y): 
	ans = (1.0/wavelength)**2 * beta(x0,y0,I0,b,kx,ky,x,y,theta)*I0*(x-x0)*(-8.0*pi**2 - 9.8696*(x-x0) * (-2.0*(x-x0)*(kx**2 * cos(theta)**2 + ky**2 * sin(theta)**2) - (kx**2 - ky**2)*(y-y0)*sin(2*theta)))
	return ans

def dy(x,y):
	ans = -9.8696*((1.0/wavelength)**2) * I0*beta(x0,y0,I0,b,kx,ky,x,y,theta) * (-2.0*(ky**2 * cos(theta)**2 + kx**2 * sin(theta)**2)*(y-y0) - (kx**2 - ky**2)*(x-x0)*sin(2.0*theta)) * (x-x0)**2
	return ans

def dx2(x,y):
	ans = (1.0/wavelength)**2 * I0*beta(x0,y0,I0,b,kx,ky,x,y,theta) * (78.9568-19.7392*(x-x0)**2 * (kx**2 * cos(theta)**2 + ky**2 * sin(theta)**2) + 39.4784*(x-x0)*(-2.0*(x-x0)*(kx**2 * cos(theta)**2 + ky**2 * sin(theta)**2) - (kx-ky)*(y-y0)*sin(2.0*theta)) + 2.4674*(x-x0)**2 * (2.0*(x-x0)*(kx**2 * cos(theta)**2 + ky**2 * sin(theta)**2)+(kx-ky)*(y-y0)*sin(2.0*theta))**2)
	return ans                                                                                                   

def dy2(x,y):
	ans = (2.0*pi/wavelength)**2 * I0*beta(x0,y0,I0,b,kx,ky,x,y,theta)*(x-x0)**2 *(-0.5*(ky**2 * cos(theta)**2 + kx**2 * sin(theta)**2) + 0.0625*(2.0*(y-y0)*(ky**2 * cos(theta)**2 + kx**2 * sin(theta)**2) + (kx-ky)*(x-x0)*sin(2*theta))**2)
	return ans 
   
def dxy(x,y):
	answer = (1.0/wavelength)**2 * I0*beta(x0,y0,I0,b,kx,ky,x,y,theta) * (x-x0)*(-9.8696*(kx-ky)*(x-x0)*sin(2.0*theta) + 19.7392*(-2.0*(y-y0)*(ky**2 * cos(theta)**2 + kx**2 * sin(theta)**2) + (-kx+ky)*(x-x0)*sin(2.0*theta)) + 2.4674*(x-x0)*(-2.0*(y-y0)*(ky**2 * cos(theta)**2 + kx**2 * sin(theta)**2) + (-kx+ky)*(x-x0)*sin(2.0*theta))*(-2.0*(x-x0)*(kx**2 * cos(theta)**2 + ky**2 * sin(theta)**2) - (kx-ky)*(y-y0)*sin(2.0*theta)))
 	return answer



xsum = 0.0
x2sum = 0.0
ysum = 0.0
y2sum = 0.0
prsum = 0.0
graphcount = 100    # number of graphs that are being created

for j in range(graphcount):
    cen = 0.0
    Ssum = S.sum()
    for y in range(imagesize):
        for x in range(imagesize):
            cen = cen + S[y,x]*x
    xcen = cen/Ssum
for k in range(graphcount):
	cen = 0.0
	Ssum = S.sum()
	for y in range(imagesize):
        	for x in range(imagesize):
            		cen = cen + S[y,x]*y
	ycen = cen/Ssum
old = np.matrix([[xcen],[ycen]])      # our initial guess for the pixel
print old

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
	while diff1 > 0.01 or diff2 > 0.01:             
        	x0 = old[0,0]
        	y0 = old[1,0]
        	for y in range(imagesize):
           		for x in range(imagesize):
                		DX += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y,theta))-1)*dx(x,y))
                		DY += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y,theta))-1)*dy(x,y))
                		DX2 += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y,theta))-1)*dx2(x,y)) - (((dx(x,y))**2)*S[y,x]/((psf(x0,y0,I0,b,kx,ky,x,y,theta))**2))
                		DXY += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y,theta))-1)*dxy(x,y)) - (dx(x,y)*dy(x,y)*S[y,x]/(psf(x0,y0,I0,b,kx,ky,x,y,theta))**2)
                		DY2 += (((S[y,x]/psf(x0,y0,I0,b,kx,ky,x,y,theta))-1)*dy2(x,y)) - (((dy(x,y))**2)*S[y,x]/((psf(x0,y0,I0,b,kx,ky,x,y,theta))**2))

        	J = np.matrix([[DX2,DXY],[DXY,DY2]])
		F = np.matrix([[DX],[DY]])
        	J2 = J.getI()
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
	prsum = (old[0,0] * old[1,0]) + prsum

xavg = xsum/graphcount
x2avg = x2sum/graphcount
yavg = ysum/graphcount
y2avg = y2sum/graphcount
avgpr = prsum/graphcount
print xavg
print yavg

standarddeviationx = sqrt(x2avg - (xavg)**2)
standarddeviationy = sqrt(y2avg - (yavg)**2)

print standarddeviationx , standarddeviationy
cov = avgpr - (xavg*yavg)
varx = standarddeviationx**2
vary = standarddeviationy**2

M1 = np.matrix([[varx,cov],[cov,vary]])
print M1
print linalg.eigvals(M1)

