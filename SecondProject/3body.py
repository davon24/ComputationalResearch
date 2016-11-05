#This program shows how to use the Velocity Verlet method to model planetary motion.
#We include earth and Jupiter

from numpy import *
import matplotlib.pyplot as plt

def acceleration(variables,params,t):
	xE=variables[0]
	yE=variables[1]
	xJ=variables[2]
	yJ=variables[3]
	MSG = params[0]
	MJG = params[1]
	MEG = params[2]
	
	rES=sqrt(xE**2+yE**2) #Earth-Sun distance
	rJS=sqrt(xJ**2+yJ**2) #Earth-Jupiter distance
	rEJ=sqrt((xJ-xE)**2+(yJ-yE)**2)
	
	#Sun's contribution to earth's acceleration
	aES = array([-xE/rES,-yE/rES])*MSG/(rES**2)
	#Sun's contribution to Jupiter's acceleration
	aJS = array([-xJ/rJS,-yJ/rJS])*MSG/(rJS**2)
        #Jupiter's contribution to earth's acceleration
	aEJ = array([-(xE-xJ)/rEJ,-(yE-yJ)/rEJ])*MJG/(rEJ**2)
        #Earth's contribution to Jupiter's acceleration
	aJE = -aEJ*MEG/MJG #From Newton's Third Law
	
        aE = aES+aEJ
        aJ = aJS+aJE
        
        return array([aE[0],aE[1],aJ[0],aJ[1]])
	
#Variables that we'll need
xE0 = 1 #This is the initial distance, in units of AU
yE0 = 0 #We start on the x axis
vEx0 = 0 #We start at aphelion or perihelion
vEy0 = 2*pi/xE0**0.5
xJ0 = 1.2*1.5**0.666666 #This is the initial distance, in units of AU
yJ0 = 0 #We start on the x axis
vJx0 = 0 #We start at aphelion or perihelion
vJy0 = 2*pi/xJ0**0.5 #Circumference over period
MSG = 4*pi**2 #Solar MG in units of AU^3/year^2
MJG = MSG/1047 #Same for Jupiter
MEG = MSG/1047
params = array([MSG, MJG, MEG],dtype='float')

dt = 0.0005 #In units of years
onehalfdtsquared = 0.5*dt**2
tmax = 20
numtimes = int(tmax/dt)

#Initialize arrays
times = linspace(0,tmax,numtimes)
#Position and velocity arrays
coordinates = zeros([numtimes,4]) #Careful with the brackets for 2D arrays!
velocities = zeros([numtimes,4])

coordinates[0] = array([xE0,yE0,xJ0,yJ0])
velocities[0]  = array([vEx0,vEy0,vJx0,vJy0])

aold = acceleration(coordinates[0],params,0)

#The main loop
for t in range(1,numtimes):
    coordinates[t] = coordinates[t-1]+velocities[t-1]*dt+onehalfdtsquared*aold
    anew = acceleration(coordinates[t],params,times[t-1])
    velocities[t] = velocities[t-1]+0.5*(anew+aold)*dt
    aold = anew #Why calculate two accelerations per iteration when I can recycle?

#Time to show results
plt.figure()
plt.plot(coordinates[:,0],coordinates[:,1])
plt.plot(coordinates[:,2],coordinates[:,3])
plt.axis('equal')
plt.show()