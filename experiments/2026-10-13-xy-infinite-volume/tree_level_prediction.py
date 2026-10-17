"""The supplied note pp2/5: (2D-1)v1(beta,epsilon)=1."""
import numpy as np
from scipy.optimize import brentq
nodes,weights=np.polynomial.legendre.leggauss(128)

def v1(beta,epsilon):
    cosine=np.cos(epsilon*nodes);w=weights*np.exp(beta*(cosine-1))
    return float(np.dot(w,cosine)/w.sum())

def endpoint(D):
    return brentq(lambda e:np.sin(e)/e-1/(2*D-1),1e-8,np.pi,xtol=1e-13)/np.pi

def beta_critical(D,q):
    target=1/(2*D-1);epsilon=np.pi*q;low=v1(0,epsilon)-target
    if abs(low)<1e-13:return 0.
    if low>0:return np.nan
    high=1.
    while v1(high,epsilon)<target:high*=2
    return brentq(lambda b:v1(b,epsilon)-target,0,high,xtol=1e-13)
