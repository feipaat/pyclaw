import scipy.optimize
from numpy import * 
from scipy import integrate
import numpy as np
import scipy
#%matplotlib inline
import matplotlib.pyplot as plt
#import mpld3
#mpld3.enable_notebook()

def test(N,tfinal,q):
    
    #Parameteres
    
    # N grid points
    # tfinal; Break time: 0.15914
    # q norm
    
    ## CLAWPACK Burgers
    from clawpack.pyclaw import examples
    claw = examples.burgers_1d.setup(mx=N)     
    claw.tfinal = tfinal
    claw.run()
    Numeric_1_Before=claw.frames[10].q[0,:]
    from clawpack.pyclaw import examples
    claw = examples.burgers_1d.setup(mx=2*N)     
    claw.tfinal = tfinal
    claw.run()
    Numeric_2_Before=claw.frames[10].q[0,:]

    ## Semianalytic Solution

    # Char.
    def func(xi):
        return (np.sin(2.0*np.pi*xi)+0.50)*(tfinal-0)+xi-x

    xx=(np.asarray(range(1,2*N,2)))/(2.0*N)
    xi = zeros(len(xx))
    sol = []
    for guess in range(N):
        for i in range(guess):
            x = xx[i]
            xi[i] = scipy.optimize.fsolve(func, xi[i-1])
        xi[N-1]=xi[0]+1
        for j in np.arange(len(xx)-2, guess-1,-1):
            x = xx[j]
            xi[j] = scipy.optimize.fsolve(func, xi[j+1]) 
        q_0_xi=np.sin(2.0*np.pi*xi)+0.5
        sol.append(np.abs(0.50-scipy.integrate.trapz(q_0_xi,xx)))

    (m,index1) = min((index2,index1) for index1,index2 in enumerate(sol))
    h1=1.0/N
    Semi_1_Before=q_0_xi

    xx=(np.asarray(range(1,4*N,2)))/(4.0*N)
    xi = zeros(len(xx))
    for i in range(guess):
        x = xx[i]
        xi[i] = scipy.optimize.fsolve(func, xi[i-1])
    xi[2*N-1]=xi[0]+1
    for j in np.arange(len(xx)-2, guess-1,-1):
        x = xx[j]
        xi[j] = scipy.optimize.fsolve(func, xi[j+1]) 
    q_0_xi=np.sin(2.0*np.pi*xi)+0.5
    
    h2=1.0/(2.0*N)
    Semi_2_Before=q_0_xi
    
    max1=np.max(np.abs(Numeric_1_Before-Semi_1_Before))
    max2=np.max(np.abs(Numeric_2_Before-Semi_2_Before))
    import math 
    p_max=math.log(max1/max2,2)
    print '\n'
    print("In the maximum norm the order of the convergence is %s" % p_max)

    # Error and convergence rate
    if q > 0:
        E_h=(h1*np.sum(np.abs(Numeric_1_Before-Semi_1_Before)**q))**(1.0/q)
        E_h_per_2=(h2*np.sum(np.abs(Numeric_2_Before-Semi_2_Before)**q))**(1.0/q)
        R_h=(E_h)/(E_h_per_2)
        import math 
        p=math.log(R_h,2)
        print '\n'
        print("In the %s norm the order of the convergence is %s" % (q, p))
    else:
        print "Error"
    if tfinal <= 1.0/(2*np.pi):
        tfinal
    elif tfinal > 1.0/(2*np.pi):
        print "After the break time"
    else:
        print "Error"
        
    fig, ax = plt.subplots()
    ax.plot(xx,q_0_xi,'og',label=r"Semianalytic")
    ax.plot(claw.frames[10].grid.p_centers[0],claw.frames[10].q[0,:],'-r', label=r"CLAWPACK")
    ax.legend(loc='lower left') # upper left corner
    ax.set_xlabel(r'Interval')
    ax.set_ylabel(r'Function value')

if __name__=="__main__":
	test(100,0.1,1)
