"""Finite-size long-distance spin quotients at measured common xi/L crossings."""
import numpy as np
from scipy.interpolate import PchipInterpolator, PPoly

def roots(q,d,kind='linear'):
    q=np.asarray(q);d=np.asarray(d)
    if not np.isfinite(d).all():return np.array([])
    if kind=='linear':
        out=list(q[d==0])
        for i in range(len(q)-1):
            if d[i]*d[i+1]<0:out.append(q[i]-d[i]*(q[i+1]-q[i])/(d[i+1]-d[i]))
        # An identically zero segment is nonunique, not a unique crossing.
        if np.any((d[:-1]==0)&(d[1:]==0)):return np.array([q[0],q[-1]])
    else:raise ValueError('Use joint_root for PCHIP')
    return np.unique(out)

def joint_root(q,Rsmall,Rlarge,kind):
    if kind=='linear':return roots(q,Rsmall-Rlarge)
    if kind!='pchip':raise ValueError('Unknown interpolation')
    if not np.isfinite(Rsmall).all() or not np.isfinite(Rlarge).all():return np.array([])
    a=PchipInterpolator(q,Rsmall);b=PchipInterpolator(q,Rlarge)
    r=PPoly(a.c-b.c,a.x).roots(extrapolate=False)
    if np.isnan(r).any():return np.array([q[0],q[-1]])
    r=np.sort(r[np.isfinite(r)&(r>=q[0])&(r<=q[-1])]);return r[np.r_[True,np.diff(r)>1e-12]] if len(r) else r

def quotient(q,small,large,kind='linear',ratio=2):
    """Columns R,logGquarter,logGhalf,logm2; one draw or central means."""
    r=joint_root(q,small[:,0],large[:,0],kind)
    if len(r)!=1:return None,len(r)
    x=r[0]
    interp=lambda y: np.array([np.interp(x,q,y[:,j]) for j in range(4)]) if kind=='linear' else PchipInterpolator(q,y,axis=0)(x)
    s=interp(small);l=interp(large)
    delta=-(l[1:]-s[1:])/(2*np.log(ratio))
    return np.r_[x,(s[0]+l[0])/2,delta],1
