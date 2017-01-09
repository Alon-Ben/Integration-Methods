#Assignment #5: Comparison of Integration Methods: Gareth Sykes

from math import sin, sqrt, cos, log10
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


#function to be integrated

def f(x):
    return (sin(sqrt(100*x)))**2


#computes the Area under the function to be integrated using the trapozoid rule

def trapIntegral(a,b,N,dx):

    x = a
    x_i = x + dx
    dA = 0.0

    for i in xrange(1,N+1):
        dA += (dx/2) * (f(x_i) + f(x))
        x = x + dx
        x_i = x_i + dx
        

    #print "Area = ", dA, " with Error = ", e, " for ", N, " slices"

    return dA

a = 0.0
b = 1.0
e = 100000

N = 1
dx = abs((b-a)/N)
In = trapIntegral(a,b,N,dx)

N = 2
dx = abs((b-a)/N)
In1 = trapIntegral(a,b,N,dx)

# initalize bounds of integral, error, and number of slices
a = 0.0
b = 1.0
N = 1
e = 100000

dx = abs((b-a)/N)

# Will stop computing the Area once the error is approx < 10-6
while e > 10e-6:
    dx = (abs(b-a)/N)
    e = abs((1.0/3)*(In1-In))
    
    plt.plot(log10(N),log10(1/e), 'ro') # x-axis: Log10(N) y-axis: Log10(1/e), Trapozoid rule plotted in red
    
    In = In1
    In1 = trapIntegral(a,b,N,dx)
    N *= 2


#ROMBERG INTEGRATION


R = np.zeros([20,20])

#Initialize R11, and R21

N = 1
dx = (abs(b-a)/N)
R11 = trapIntegral(a,b,N,dx)

N = 2
dx = (abs(b-a)/N)
R21 = trapIntegral(a,b,N,dx)

R[0,0] = R11
R[1,0] = R21

e = 10000 #initialize error, and indicies i, m and NumOfElements
i = 1
m = 0
NumOfElements = 1

#Starts printing the table
print R[0,0]

while e > 10e-6:
    
    for r in range(NumOfElements):
        R[i,m+1] = R[i,m] + (1.0/(4.0**(m+1)-1))*(R[i,m] - R[i-1,m])
        m += 1
        if i-1 == m :
            print R[i,m]
        else:
            print R[i,m],

    #Set up variables for the next while loop    
    m = 0
    i += 1
    NumOfElements += 1
    N *= 2
    dx = (abs(b-a)/N)
    R[i,m] = trapIntegral(a,b,N,dx)
    
    e = abs((1.0/(4.0**(NumOfElements-1)-1))*(R[NumOfElements-1,NumOfElements-1] - R[(NumOfElements-1)-1,NumOfElements-1]))
    plt.plot(log10(N/2), log10(1/e), 'bo') #N/2 since I already changed N to the next value above
    # x-axis: Log10(N) y-axis: Log10(1/e), Remberg rule plotted in blue


#SIMPSONS RULE
def simpsons(a,b,N):
    
    dx = abs((b-a)/N)

    OddSum = 0.0
    EvenSum = 0.0

    for k in range(1,N,2):
        OddSum += f(a+k*dx)

    for k in range(2,N,2):
        EvenSum += f(a+k*dx)

    I = (1.0/3.0)*dx*(f(a) + f(b)+4*OddSum + 2*EvenSum)

    #print I, " Simspons Rule for ", N, " slices"

    return I

N = 2 #initialize N to 2
In = simpsons(a,b,2)
In1 = simpsons(a,b,4)
SimpError = 100000

while (SimpError > 10e-6):
    simpsons(a,b,N)
    
    SimpError = abs((1.0/15)*(In1-In))
    #print "Error on Simp = ", SimpError, " N = ", N
    plt.plot(log10(N), log10(1/SimpError), 'go') # x-axis: Log10(N) y-axis: Log10(1/SimpError), Simpsons rule plotted in green
    
    In = In1
    In1 = simpsons(a,b,N)
    N *= 2




#PLOTING OF ERRORS ( inbedded in while loops for each method ):


plt.savefig('figure.png')

#The Required accuracy is log(1/10e-6) or 5 --> simpsons rule and Romberg Integration reach this accuracy at about the same time ( as seen on the graph ) with the same number of slices, however, from the data that I have, simpsons rule is slightly faster

# it's also interesting to note that different methods have different levels of accuracy relative to the other methods on specific intervals ( ie simpsons has the greatest accuracy at the beggining but is quickly overtaken by Romberg Integration then Simpsons rule surpasses Romberg Integration again at Log(N) ~ 5.
