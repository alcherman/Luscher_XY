import itertools, math
import numpy as np
from supervillain import Lattice,XY
from supervillain.generator.xy import ClusterUpdate

def distance(a,b):
    return abs(math.atan2(math.sin(a-b),math.cos(a-b)))

def check(D,N,kappa,epsilon,alpha):
    shape=(N,)*D
    coords=list(itertools.product(range(N),repeat=D))
    coord_index={c:i for i,c in enumerate(coords)}
    edges=[]
    for mu in range(D):
        for i,c in enumerate(coords):
            d=list(c);d[mu]=(d[mu]+1)%N
            edges.append((mu,i,coord_index[tuple(d)]))
    V=len(coords)
    # Branch-cut crossing and mixed reflection-projection signs.
    base=(np.linspace(-.5,.5,V)+5.9)%(2*math.pi)
    rb=np.angle(-np.exp(2j*alpha)*np.exp(-1j*base))%(2*math.pi)
    bits=np.array(list(itertools.product([0,1],repeat=V)))
    fields=np.where(bits,rb,base)
    E=[]
    for f in fields:
        E.append(-kappa*sum(math.cos(f[i]-f[j]) for mu,i,j in edges)
                 if all(distance(f[i],f[j])<epsilon for mu,i,j in edges) else math.inf)
    allowed=np.flatnonzero(np.isfinite(E));indexes={v:i for i,v in enumerate(allowed)}
    P=np.zeros((len(allowed),len(allowed)))
    action=XY(Lattice(D,N),kappa,epsilon);update=ClusterUpdate(action)
    graphs=moves=0
    for s in allowed:
        f=fields[s];r=np.angle(-np.exp(2j*alpha)*np.exp(-1j*f))%(2*math.pi)
        ps=np.array([1. if distance(r[i],f[j])>=epsilon else
            -math.expm1(-max(0.,kappa*(math.cos(f[i]-f[j])-math.cos(r[i]-f[j]))))
            for mu,i,j in edges])
        variable=np.flatnonzero((ps>0)&(ps<1))
        for choices in itertools.product([False,True],repeat=len(variable)):
            on=ps==1;on[variable]=choices
            q=np.prod(np.where(on,ps,1-ps))
            tape=np.minimum(np.where(on,ps/2,(1+ps)/2),np.nextafter(1.,0.)).reshape(D,V)
            graphs+=1
            for seed in range(V):
                result=np.asarray(update.cluster(action.canonical(f.reshape(shape)),seed,alpha,tape)).ravel()
                b=(np.abs(result-rb)<np.abs(result-base)).astype(int)
                t=int(b @ (2**np.arange(V-1,-1,-1)))
                assert t in indexes
                assert np.max(np.abs(result-fields[t]))<1e-14
                P[indexes[s],indexes[t]]+=q/V;moves+=1
    W=np.exp(-np.array(E)[allowed]);W/=W.sum()
    err=np.max(np.abs(W[:,None]*P-(W[:,None]*P).T))
    assert np.max(np.abs(P.sum(1)-1))<5e-14
    assert err<5e-14
    print((D,N,kappa,epsilon,alpha),'allowed',len(allowed),'graphs',graphs,'moves',moves,'max_flux_error',err)
for args in [(1,3,.8,1.8,.73),(1,3,.9,math.pi,1.8),(2,2,.7,1.9,.13),(2,2,0,1.9,.13),(2,2,.5,math.pi,2.6)]:
    check(*args)
