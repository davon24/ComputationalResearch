from numpy import *
import matplotlib.pyplot as plt




def derivatives(vars, params, t):
    #params is an array of parameters. 
    #t is the time
  
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
      
#Now we have to use the function "derivatives"
  
#first initialize a bunch of needed variables
#We start with things relating to time

dt = 50.0/100.0
w = 2.0*pi
R = 1.0
msg = 4.0 * pi**2
meg = msg/100.0
params = array([w,msg,R,meg])
#Now set up the variables that we'll need
numtimes= 3*10**4
times=zeros(numtimes)
xs0 = 0.0
ys0 = 13.0
vxs0 = sqrt(msg/ys0)
vys0 = 0.0
  
xs_old = xs0
ys_old = ys0
vsx_old = vxs0
vsy_old = vys0

#Now we do the Runge-Kutta Method
i=0
j=0


xs = [xs0]
ys = [ys0]
xs_poin = zeros(numtimes)
vsx_poin = zeros(numtimes)
vsy_poin = zeros(numtimes)
outfile1 = open('solarinterpolated6.txt' , 'w') #opens file where we will store our values
outfile1.write('The interpolated values for the x position and x velocity are in the first and second column respectively.' '\n' 'There are %d points. The spaceship starts at (%f,%f) and the Earth starts at (1,0). The time step is %f' '\n''\n' % (numtimes, xs0, ys0, dt))
outfile2 = open('solarnoninterpolated6.txt' , 'w')
outfile2.write('The interpolated values for the x position and x velocity are in the first and second column respectively.' '\n' 'There are %d points. The spaceship starts at (%f,%f) and the Earth starts at (1,0). The time step is %f' '\n''\n' % (numtimes, xs0, ys0, dt))
outfile3 = open('shiptrajectory6.txt', 'w')
outfile3.write('The interpolated values for the x position and x velocity are in the first and second column respectively.' '\n' 'There are %d points. The spaceship starts at (%f,%f) and the Earth starts at (1,0). The time step is %f' '\n''\n' % (numtimes, xs0, ys0, dt))
while j<numtimes: #runs while loops that calculates and stores the points to be used in the poincare section
	

	x = array([xs_old,ys_old,vsx_old,vsy_old])

	#The next four lines calculate derivatives at different times, and multiply the derivatives by the time step, to get the changes in the variables
	k1 = dt*derivatives(x,params,dt*(i))
	k2 = dt*derivatives(x+k1/2,params,dt*(i)+dt/2)
	k3 = dt*derivatives(x+k2/2,params, dt*(i)+dt/2)
	k4 = dt*derivatives(x+k3,params, dt*(i)+dt)
	#Now that we've done the RK method and solved the ODE with the old time and data at the previous step, we up the time and step and store for the new time
	x_new = x+k1/6+k2/3+k3/3+k4/6 #This is the Runge-Kutta estimate of the new values of theta and d(theta)/dt
	[xs_new,ys_new,vsx_new,vsy_new]=x_new
	
	i+=1
	xs.append(xs_new)
	ys.append(ys_new)
	
	outfile3.write('%4f  %4f\n' %(xs[i],ys[i]))
	if ys_new <= 0 and ys_old > 0:
		
		s = ys_new/(ys_new - ys_old)
		xs_poin[j] = s*xs_old + (1-s)*xs_new
		vsx_poin[j] = s*vsx_old + (1-s)*vsx_new
		vsy_poin[j] = s*vsy_old +(1-s)*vsy_new
		outfile1.write('%4f  %4f  %4f\n' % (xs_poin[j], vsx_poin[j], vsy_poin[j]))
		outfile2.write('%4f  %4f  %4f\n' %(xs_new, vsx_new, vsy_new))
		j += 1
		if j/100.0 == int(j/100.0):
			print('We have found %d out of %d points for the Poincare section.' % (j, numtimes))
	

	xs_old = xs_new
	ys_old = ys_new
	vsx_old = vsx_new
	vsy_old = vsy_new
	if i/(j+1) > 170:
		break
#outfile1.close()
#outfile2.close()
  
 

#next we check our results

plt.figure(1)
plt.plot(xs,ys, 'bo')
plt.show()
