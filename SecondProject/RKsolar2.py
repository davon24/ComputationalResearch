from numpy import *
import matplotlib.pyplot as plt




def derivatives(vars, params, t):
    #params is an array of parameters. 
    #t is the time

	xs = vars[0]
	ys = vars[1]
	vsx = vars[2]
	vsy = vars[3]
	
	w = params[0]    
	msg = params[1]
	R = params[2]
	meg = params[3]

	xe = R * cos(w*t) #The x coordinate of the earth's position
	ye = R * sin(w*t) #The y coordinate of the earth's position
	Rss3 = (xs**2 + ys**2)**1.5
	Fxss = -(xs*msg)/Rss3 #Force between the sun and the satellite in the x direction
	Fyss = -(ys*msg)/Rss3 #Force between the sun and the satellite in the y direction
	
	Rse3 = ((xs - xe)**2 + (ys - ye)**2)**1.5
	Fxse = -(meg)*(xs-xe)/Rse3 #Force between the earth and the satellite in the x direction
	Fyse = -(meg)*(ys-ye)/Rse3 #Force between the earth and the satellite in the y direction
	
	xdot = vsx
	ydot = vsy
	
	vxdot = Fxss + Fxse #Acceleration of the satellite in the x direction
	vydot = Fyss + Fyse #Acceleration of the satellite in the y direction
	return array([xdot, ydot, vxdot, vydot])
      
#Now we have to use the function "derivatives"
  
#first initialize a bunch of needed variables
#We start with things relating to time

dt = 1.0/200.0
w = 2.0*pi
R = 5.2
msg = 4.0 * pi**2
meg = msg/10.0
params = array([w,msg,R,meg])
#Now set up the variables that we'll need
numtimes= 1000
times=zeros(numtimes)
xs0 = 0.0
ys0 = 3.276
vxs0 = 0
vys0 = 3.471
  
xs_old = xs0
ys_old = ys0
vsx_old = vxs0
vsy_old = vys0

#Now we do the Runge-Kutta Method
i=0
j=0


xs = [xs0]
ys = [ys0]


while i<numtimes*300: 
	

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
	
	
	xs_old = xs_new
	ys_old = ys_new
	vsx_old = vsx_new
	vsy_old = vsy_new
     

  
 

#next we check our results

plt.figure(1)
plt.plot(xs,ys, 'bo')
plt.show()
