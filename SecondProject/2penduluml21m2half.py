#This program simulates a double pendulum for different initial conditions

#This version implements adaptive step size

from numpy import *
import matplotlib.pyplot as plt

def derivatives(vars,params,t):
    theta1 = vars[0]
    theta2=vars[1]
    p1=vars[2]
    p2=vars[3]
    m2=params[0] #We assume that m1=1
    l2=params[1] #We assume that l1=1
    g = params[2]	
	
    mess = l2*(1+m2*(sin(theta1-theta2))**2)
    #We need this, so we pre-compute it

    thetadot1 = (l2*p1-p2*cos(theta1-theta2))/mess
    thetadot2 = ((1+m2)*p2-l2*m2*p1*cos(theta1-theta2))/(mess*m2*l2)
	
    c1 = (p1*p2*sin(theta1-theta2))/(l2*(1+m2*(sin(theta1-theta2))**2))
    c2 = (m2*((l2)**2)*(p1**2)) + (p2**2)*(1.0+m2) - 2*l2*m2*p1*p2*cos(theta1-theta2)
    c2 = c2/(2.0*(l2**2)*(1.0+m2*sin(theta1-theta2)**2)**2)
    c2 = c2*sin(2.0*(theta1-theta2))

    pdot1 = -(1+m2)*g*sin(theta1) - c1 + c2
    pdot2 = -m2*g*sin(theta2) + c1 - c2
	
    return array([thetadot1,thetadot2,pdot1,pdot2])

def Hamiltonian(theta1,theta2,p1,p2,params):
    m2=params[0]
    l2=params[1]
    g=params[2]
    
    KE = (m2*(l2*p1)**2 + (1+m2)*p2**2 - 2*m2*l2*p1*p2*cos(theta1-theta2))/(2*m2*(l2**2)*(1+m2*(sin(theta1-theta2))**2))
    U = g*((1+m2)*(1-cos(theta1))+m2*l2*(1-cos(theta2)))
    
    return KE+U

#Variables that we'll need
p10=0.0 #We assume that each pendulum starts from rest
p20=0.0
m2=0.5 #We assume that m1 = 1, the only variable is m2
l2=1.0 #We assume that l1 = 1, the only variable is m2
g=4*pi**2 #We picked this value so the period of small oscillations of a single pendulum would be 1
params=array([m2,l2, g])
log10 = log(10) #We'll be using this over and over, to convert times to a logarithmic scale


#Compute the minimum energy needed to flip
energyflip = Hamiltonian(0,pi,0,0, params)

maxtime=5*10**2
numangles1 = 500
numangles2 = 2*numangles1-1
thetamax = pi-0.01
dt0 = 0.25 #We'll adapt the step size, but this is what we'll start with.  1/400 of a cycle of small oscillation is good.
dtmax = dt0*10 #We'll never let the step size get above this
errortarget = 10**(-8) #We want an error no larger than a specified fraction of the energy


#Start the loops
outfile=open('logtflipl21m2half.dat', 'w')
for i in range(0,numangles1):
    theta10 = thetamax*i/(numangles1-1) #We vary theta10 between 0 and thetamax
    print 'i = ', i, ', starting value of theta1 is ', theta10
   
    for j in range(0,numangles2):
        theta20 = 2*thetamax*j/(numangles2-1)-thetamax #We vary theta20 between -thetamax and +thetamax
        
        theta2=theta20
        p1=p10
        p2=p20
        theta1=theta10
 
        t=0 #Start at t=0
        dt = dt0 #We'll increase or decrease this as needed
        #We will speed things up 50% by tossing out data points for which the energy is too low to flip
        energyi = Hamiltonian(theta1,theta2,p1,p2,params)
        keepgoing=True
        if energyi<energyflip:
            t=1.1*maxtime #This way we'll never enter the loop below, and we clearly mark these cases in the data
        
        while(t<maxtime and abs(theta2)<pi and keepgoing):
            #Runge-Kutta 4th order for a half-step
            k1 = 0.5*dt*derivatives(array([theta1,theta2,p1,p2]),params,t)
            k2 = 0.5*dt*derivatives(array([theta1,theta2,p1,p2])+k1/2,params,t+dt/4)
	    k3 = 0.5*dt*derivatives(array([theta1,theta2,p1,p2])+k2/2,params,t+dt/4)
	    k4 = 0.5*dt*derivatives(array([theta1,theta2,p1,p2])+k3,params,t+dt/2)
	    [theta1a,theta2a,p1a,p2a]=array([theta1,theta2,p1,p2])+k1/6+k2/3+k3/3+k4/6

            #Check the angles half-way through when we do the double-step
	    if abs(theta2a)>pi:
	        keepgoing = False
	        #It's  possible that the first half-step will send it past pi, then the next half-step will send it back.
	        #Using a full step would never find that. This way we spot it.                        
            
            #Next half-step
            k1 = 0.5*dt*derivatives(array([theta1a,theta2a,p1a,p2a]),params,t+dt/2)
            k2 = 0.5*dt*derivatives(array([theta1a,theta2a,p1a,p2a])+k1/2,params,t+3*dt/4)
	    k3 = 0.5*dt*derivatives(array([theta1a,theta2a,p1a,p2a])+k2/2,params,t+3*dt/4)
	    k4 = 0.5*dt*derivatives(array([theta1a,theta2a,p1a,p2a])+k3,params,t+dt)
	    [theta1b,theta2b,p1b,p2b]=array([theta1a,theta2a,p1a,p2a])+k1/6+k2/3+k3/3+k4/6            
            
            #Now do a full step
            k1 = dt*derivatives(array([theta1,theta2,p1,p2]),params,t)
            k2 = dt*derivatives(array([theta1,theta2,p1,p2])+k1/2,params,t+dt/2)
	    k3 = dt*derivatives(array([theta1,theta2,p1,p2])+k2/2,params,t+dt/2)
	    k4 = dt*derivatives(array([theta1,theta2,p1,p2])+k3,params,t+dt)
	    [theta1c,theta2c,p1c,p2c]=array([theta1,theta2,p1,p2])+k1/6+k2/3+k3/3+k4/6
	    
            t+=dt
            theta1=theta1b
            theta2=theta2b
            p1=p1b
            p2=p2b
            #Let's see if energy is conserved
            energy = Hamiltonian(theta1,theta2,p1,p2,params) #Energy from the pair of half-steps
       	    
       	    #Now we compute our error by comparing the two approaches, and pick a new timestep
       	    energyfull = Hamiltonian(theta1c,theta2c,p1c,p2c,params)
       	    error  = abs((energyfull-energy)/energy)
            if error>0:
                dt = dt*(errortarget/error)**0.2
            if dt>dtmax:
                dt = dtmax

       
       	#Let's see if energy was conserved
       	energy = Hamiltonian(theta1,theta2,p1,p2,params)
       	#We recompute because if it exited the loop for insufficient energy we'll be comparing energyi with old energy
       	#print 'iteration', i*numangles2+j+1, 'of ', numangles1*numangles2, ', theta10=', theta10, ', theta20=', theta20, ', fliptime=', t, ', energy loss = ', (1-energy/energyi)
        #Output the latest calculations to a file
        outfile.write('%14.8f' % (log(t)/log10))
    #After we've finished varying theta20 we add a line break to the data file
    outfile.write('\n')

#Need to output a file
outfile.close()
