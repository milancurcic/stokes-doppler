#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('k_long',help='swell wavenumber')
parser.add_argument('a_long',help='swell amplitude')
parser.add_argument('k_short',help='windsea wavenumber')
args = parser.parse_args()

k_l = float(args.k_long)
a = float(args.a_long)
k0 = float(args.k_short)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def u_l(x,z,t,a,k,omega):
    """
    Returns the long wave horizontal velocity at input level z.
    """
    return a*omega*np.sin(k*x-omega*t)*np.exp(k*z)

def stokes_drift(x,a,k,omega):
    """
    Returns the surface Stokes drift.
    """
    return np.mean(u_l(x,eta(x,0,a,k,omega),0,a,k,omega))

def eta(x,t,a,k,omega):
    """
    Returns the long wave surface elevation.
    """
    return a*np.sin(k*x-omega*t)

def diff(x,periodic=False):
    """
    1-dimensional 2-nd order centered differencing.
    """
    dx = np.zeros(x.size)
    dx[1:-1] = 0.5*(x[2:]-x[0:-2])
    if periodic:
        dx[0] = 0.5*(x[1]-x[-1])
        dx[-1] = 0.5*(x[0]-x[-2])
    else:
        dx[0] = x[1]-x[0]
        dx[-1] = x[-1]-x[-2]
    return dx

def omega(g,k,u,sigma=0.07,rho=1e3):
    """
    Returns Doppler-shifted angular frequency.
    """
    return np.sqrt(g*k + sigma/rho*k**3) + k*u

def integrand(k,t):
    return -diff(omega(g,k,u_l(x,eta(x,t,a,k_l,np.sqrt(g*k_l)),t,a,k_l,np.sqrt(g*k_l))),periodic=True)/diff(x)

def crest_position(x,t,a,k,omega):
    return np.array([x[np.argmax(eta(x,t[n],a,k,omega))] for n in range(t.size)])

def rk4( f, x0, t ):
    """
    Fourth-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.

    AUTHOR:
        Jonathan Senning <jonathan.senning@gordon.edu>
        Gordon College
        Based Octave functions written in the spring of 1999
        Python version: March 2008, October 2008

    USAGE:
        x = rk4(f, x0, t)

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a NumPy array.  In this
                case f must return a NumPy array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or NumPy array
                if a system of equations is being solved.
        t     - list or NumPy array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - NumPy array containing solution values corresponding to each
                entry in t array.  If a system is being solved, x will be
                an array of arrays.
    """
    n = len( t )
    x = np.array( [ x0 ] * n )
    for i in xrange( n - 1 ):
        h = t[i+1] - t[i]
        k1 = h * f( x[i], t[i] )
        k2 = h * f( x[i] + 0.5 * k1, t[i] + 0.5 * h )
        k3 = h * f( x[i] + 0.5 * k2, t[i] + 0.5 * h )
        k4 = h * f( x[i] + k3, t[i+1] )
        x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0

    return x


phi = np.linspace(0,2*np.pi,1e2,endpoint=False) # phase space [rad]
g = 9.8 # gravitational acceleration [m s^-2]
x = phi/k_l # physical space [m]

dt = 1e-4
duration = 10
time = np.arange(0,duration+dt,dt)
k = np.ones(x.size)*k0

# integrate short wavenumber in time
knew = rk4(integrand,k,time)

print k_l,a,np.min(knew),np.max(knew),np.mean(knew)

# calculate the position of swell crest in physical space and recast to phase space
crest = crest_position(x,time,a,k_l,np.sqrt(g*k_l))*k_l

kdiff = knew/k*100-100
kdiff[kdiff >  29.999] =  29.999
kdiff[kdiff < -29.999] = -29.999

# Evaluate short wave group velocity and propagation line
cg = 0.5*np.sqrt(g/k0)
crest_short = 0.5*np.pi+cg*time*k_l

# Evaluate long wave Stokes drift and propagation line
stokes_long = stokes_drift(x,a,k_l,np.sqrt(g*k_l))
stokes_line = 0.5*np.pi+stokes_long*time*k_l

cg_plus_stokes_line = 0.5*np.pi+(cg+stokes_long)*time*k_l

# Plot and save to file
fig = plt.figure(figsize=(8,7))
ax = fig.add_subplot(111,xlim=(0,2*np.pi),ylim=(0,duration))
plt.contourf(phi,time,kdiff,np.arange(-30,31,1),cmap=cm.bwr)
plt.colorbar(ticks=range(-30,35,5))
plt.plot(crest,time,'k.',ms=1)
plt.plot(crest_short,time,'k--',lw=2)
plt.plot(stokes_line,time,'c--',lw=2)
plt.plot(cg_plus_stokes_line,time,'m--',lw=2)
plt.xlabel('Phase [rad]',fontsize=16)
plt.ylabel('Time [s]',fontsize=16)
plt.xticks(np.arange(0,2.5*np.pi,0.5*np.pi))
ax.set_xticklabels([r'$0$',r'$\pi/2$',r'$\pi$',r'$3\pi/2$',r'$2\pi$'])
plt.title(r'$k_{short}$ change [%], $k_{long} = $'+'%3.1f'%k_l+', $a_{long} = $'+'%5.3f'%a+', $k_{short} = $'+'%3i'% int(k0),fontsize=16)
plt.grid(True)
plt.savefig('kshort_kl='+'%3.1f'%k_l+'_a='+'%5.3f'%a+'_ks='+'%3.3i'% int(k0)+'.png',dpi=100)
plt.close(fig)
