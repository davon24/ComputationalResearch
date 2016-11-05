from numpy import *
import matplotlib.pyplot as plt
  
  
def derivatives(vars, params, t):
    #This is to solve a pendulum problem, where vars[0] is theta and vars[1] is the time derivative of theta.
    #Most of the physics is contained in this function.
    #params is an array of parameters.  In this case, it's just g.
    #t is the time
      
    theta = vars[0]
    omega=vars[1]
    goverL=params[3]    
    q = params[0]
    w = params[1]
    F = params[2]
    omegadot = -goverL*sin(theta) - q*omega + F*sin(w*t) #q =friction F = driving force     w = driving frequency
    return array([omega,omegadot])
      
#Now we have to use the function "derivatives"
  
#first initialize a bunch of needed variables
#We start with things relating to time
ti = 0
dt = 0.5/100.0
theta_i = 0.2
goverL=1 #We pick a value so that the frequency of small oscillations is omega=2pi and the period is 1
thetadot_i = 0
q = 0.5
w = 2.0/3.0
F = 1.2
  
params = array([q,w,F,goverL])
#Now set up the variables that we'll need
numtimes= .25 * 10**5
times=zeros(numtimes)
  
  
theta_old = theta_i
omega_old = thetadot_i
 
theta_poin = zeros(numtimes)
omega_poin = zeros(numtimes)
 
#Now we do the Runge-Kutta Method
i=0
j=0
wdt = w*dt

def first_sine(i):
	ans = sin(wdt*i)
	return ans

def second_sine(i):
	ans = sin(wdt*(i-1))
	return ans
	
outfile1 = open('idata05.txt' , 'w') #opens file where we will store our values
outfile1.write('The interpolated values for theta and omega are in the first and second column respectively.' '\n' 'There are %d points, F = %f2 , frequency = %f, q = %f2 and dt = %f' '\n''\n' % (numtimes, F, w, q, dt))
outfile2 = open('nidata05.txt' , 'w')
outfile2.write('The non-interpolated values for theta and omega are in the first and second column respectively.' '\n' 'There are %d points, F = %f , frequency = %f, q = %f and dt = %f' '\n''\n' % (numtimes, F, w, q, dt))

while j<numtimes: #runs while loops that calculates and stores the points to be used in the poincare section
	
	if theta_old > pi:
		theta_old -= 2.0*pi
	elif theta_old < -pi:
		theta_old += 2.0*pi
	x = array([theta_old,omega_old])
	#The next four lines calculate derivatives at different times, and multiply the derivatives by the time step, to get the changes in the variables
	k1 = dt*derivatives(x,params,dt*(i))
	k2 = dt*derivatives(x+k1/2,params,dt*(i)+dt/2)
	k3 = dt*derivatives(x+k2/2,params, dt*(i)+dt/2)
	k4 = dt*derivatives(x+k3,params, dt*(i)+dt)
	#Now that we've done the RK method and solved the ODE with the old time and data at the previous step, we up the time and step and store for the new time
	x_new = x+k1/6+k2/3+k3/3+k4/6 #This is the Runge-Kutta estimate of the new values of theta and d(theta)/dt
 
	i+=1
	
	
	if (first_sine(i) >= 0) and (second_sine(i) <= 0) and (i*dt > 500):
		s = second_sine(i)/(second_sine(i) - first_sine(i))
		theta_poin[j] = (1-s)*theta_old + s*x_new[0]
		omega_poin[j] = (1-s)*omega_old + s*x_new[1]
		outfile1.write('%4f  %4f\n' % (theta_poin[j], omega_poin[j]))
		outfile2.write('%4f  %4f\n' %(x_new[0], x_new[1]))
		j += 1
		if j/100.0 == int(j/100.0):
			print('We have found %d out of %d points for the Poincare section.' % (j, numtimes))
	

	theta_old = x_new[0]
	omega_old = x_new[1]
     
outfile1.close()
outfile2.close()
  
 
 
#next we check our results
 
plt.figure(1)
plt.plot(theta_poin,omega_poin, 'bo')
plt.show()
