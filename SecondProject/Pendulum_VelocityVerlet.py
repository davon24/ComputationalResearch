#This program simulates a pendulum with Heun's method

from numpy import *
import matplotlib.pyplot as plt

def acceleration(coordinates,params,t):
    g=params[0]
    theta = coordinates[0]
    
    omegadot = -g*sin(theta)
    	
    return omegadot

def energy(coordinates,velocities,params,t):
    g=params[0]
    theta = coordinates[0]
    omega = velocities[0]
    
    KE = 0.5*omega**2
    U = g*(1-cos(theta))
    
    return KE+U

#Variables that we'll need
omega0=0.0 #We assume that each pendulum starts from rest
theta0 = 1.57

g=4*pi**2 #We picked this value so the period of small oscillations of a single pendulum would be 1
params=array([g])


numtimes = 2*10**3
dt = 0.0125
onehalfdtsquared = 0.5*dt**2

#Create arrays
times = linspace(0, (numtimes-1)*dt, numtimes)
energies = zeros(numtimes)
coordinates = zeros([numtimes,1])
velocities = zeros([numtimes,1])

coordinates[0] = theta0
velocities[0] = omega0
aold = acceleration(coordinates[0],params,0)
energies[0] = energy(coordinates[0],velocities[0], params,0)

#The main loop
for t in range(1,numtimes):
    coordinates[t] = coordinates[t-1]+velocities[t-1]*dt+onehalfdtsquared*aold
    anew = acceleration(coordinates[t],params,times[t-1])
    velocities[t] = velocities[t-1]+0.5*(anew+aold)*dt
    aold = anew #Why calculate two accelerations per iteration when I can recycle?
    energies[t] = energy(coordinates[t],velocities[t], params,times[t])

plt.figure()
plt.subplot(211)
plt.plot(times,coordinates[:,0],label=r'$\theta$ (Velocity Verlet)')  #Make a plot
plt.xlabel("Time")
plt.ylabel("Angle")
plt.legend(loc = 'upper right')
plt.subplot(212)
plt.plot(times,energies,label='Energy (Velocity Verlet)')
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend(loc = 'upper right')
plt.show()