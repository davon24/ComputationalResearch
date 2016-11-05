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
 
dt = 0.025/100.0
m2 = 1.0 #Assume m1 = 1
l2 = 1.0 #Assume l2 = 1
g = 4.0*pi**2
params = array([m2,l2,g])
#Now set up the variables that we'll need
numtimes= 10000
times=zeros(numtimes)
alpha0 = 90.0*(pi/180.0)
beta0 = 0.0*(pi/180.0)
L_alpha0 = 0.0
L_beta0 = 0.0
 
alpha_old = alpha0
beta_old = beta0
L_alpha_old = L_alpha0
L_beta_old = L_beta0
 
#Now we do the Runge-Kutta Method
i=0
j=0
 
'''def alpha_dot(beta,L_alpha,L_beta):
    alphadot = 2.0*(L_alpha - (1.0 + cos(beta))*L_beta)/(3.0 - cos(2.0*beta))
    return alphadot
 
def Energy(alpha, LA, beta, LB):
    energy = -2*cos(alpha) - cos(alpha + beta) + (LA**2 - 2.0*(1.0+cos(beta))*LA*LB + (3.0+2*cos(beta))*LB**2)/(3-cos(2*beta))
    return energy
 
Energytotal = Energy(alpha0, L_alpha0, beta0, L_beta0)'''
 
alpha = [alpha0]
beta = [beta0]
beta_poin = zeros(numtimes)
L_alpha_poin = zeros(numtimes)
L_beta_poin = zeros(numtimes)
 
 
outfile1 = open('00025interp.txt' , 'w') #opens file where we will store our values
outfile1.write('The interpolated values for the x and y position and value of theta 2 are in the first second and third column respectively.' '\n' 'There are %d points. The pendulum starts at %f for the first angle and %f for the second angle. The time step is %f' '\n''\n' % (numtimes, alpha_old, beta_old, dt))
outfile2 = open('00025noninterp.txt' , 'w')
outfile2.write('The interpolated values for the x and y position and value of theta 2 are in the first second and third column respectively.' '\n' 'There are %d points. The pendulum starts at %f for the first angle and %f for the second angle. The time step is %f' '\n''\n' % (numtimes, alpha_old, beta_old, dt))
 
while j<numtimes: #runs while loops that calculates and stores the points to be used in the poincare section
     
    x = array([alpha_old,beta_old,L_alpha_old,L_beta_old])
    #The next four lines calculate derivatives at different times, and multiply the derivatives by the time step, to get the changes in the variables
    k1 = dt*derivatives(x,params,dt*(i))
    k2 = dt*derivatives(x+k1/2,params,dt*(i)+dt/2)
    k3 = dt*derivatives(x+k2/2,params, dt*(i)+dt/2)
    k4 = dt*derivatives(x+k3,params, dt*(i)+dt)
    #Now that we've done the RK method and solved the ODE with the old time and data at the previous step, we up the time and step and store for the new time
    x_new = x+k1/6+k2/3+k3/3+k4/6 #This is the Runge-Kutta estimate of the new values of theta and d(theta)/dt
    [alpha_new,beta_new,L_alpha_new,L_beta_new]=x_new
    i+=1
 
    '''alpha_dot_new = alpha_dot(beta_new,L_alpha_new,L_beta_new)'''
 
     
     
    '''alpha.append(alpha_new)
    beta.append(beta_new)'''
 
 
 
    if alpha_new >= 0 and alpha_old <= 0:
         
        s = alpha_old/(alpha_old - alpha_new)
        beta_poin[j] = s*beta_new + (1-s)*beta_old
        L_alpha_poin[j] = s*L_alpha_new + (1-s)*L_alpha_old
        L_beta_poin[j] = s*L_beta_new +(1-s)*L_beta_old
        outfile1.write('%4f  %4f  %4f\n' % (beta_poin[j], L_alpha_poin[j], L_beta_poin[j]))
        outfile2.write('%4f  %4f  %4f\n' %(beta_new, L_alpha_new, L_beta_new))
        j += 1
 
        if j/100.0 == int(j/100.0):
            print('We have found %d out of %d points for the Poincare section.' % (j, numtimes))
     
    '''if alpha_old >= 0 and alpha_new <= 0:
         
        s = alpha_new/(alpha_new - alpha_old)
        beta_poin[j] = s*beta_old + (1-s)*beta_new
        L_alpha_poin[j] = s*L_alpha_old + (1-s)*L_alpha_new
        L_beta_poin[j] = s*L_beta_old +(1-s)*L_beta_new
        outfile1.write('%4f  %4f  %4f\n' % (beta_poin[j], L_alpha_poin[j], L_beta_poin[j]))
        outfile2.write('%4f  %4f  %4f\n' %(beta_new, L_alpha_new, L_beta_new))
        j += 1
 
        if j/100.0 == int(j/100.0):
            print('We have found %d out of %d points for the Poincare section.' % (j, numtimes))'''
 
    alpha_old = alpha_new
    beta_old = beta_new
    L_alpha_old = L_alpha_new
    L_beta_old = L_beta_new
    #if i/(j+1) > 170:
        #break
 
outfile1.close()
outfile2.close()
   
 
#next we check our results
'''plt.figure(1)
plt.plot(alpha, beta, 'bo')
plt.show()'''

