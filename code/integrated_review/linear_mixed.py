"""REML with nonnegative variance components and no imposed upper ceiling.
Derived from the independently checked working-likelihood fitter; kept local so
changes to this review do not silently alter earlier audit implementations.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize

def fit(y,X,S,K,start=None):
    n,p=X.shape
    if np.linalg.matrix_rank(X)!=p or n<=p:raise ValueError('Rank or residual df failure')
    Ks=np.array(K)
    def core(q,ret=False):
        C=S+np.einsum('k,kij->ij',q,Ks)+np.eye(n)*1e-9
        cf=cho_factor(C,lower=True,check_finite=False)
        inv=cho_solve(cf,np.eye(n),check_finite=False);ix=inv@X;bc=np.linalg.inv(X.T@ix);b=bc@(ix.T@y);P=inv-ix@bc@ix.T;py=P@y
        obj=2*np.log(np.diag(cf[0])).sum()+2*np.log(np.diag(np.linalg.cholesky(X.T@ix))).sum()+y@py
        grad=np.einsum('ij,kji->k',P,Ks)-np.einsum('i,kij,j->k',py,Ks,py)
        if ret:return dict(beta=b,bc=bc,C=C,P=P,inv=inv,obj=float(obj),q=q,gradient=grad)
        return obj,grad
    bounds=[(0,None) if np.any(k) else (0,0) for k in K]
    starts=[np.full(len(K),v) for v in [0.,.1,1.,10.]]
    if start is not None:starts.append(np.asarray(start))
    opts=[]
    for x in starts:
        opts.append(minimize(core,x,jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':400,'ftol':1e-10,'gtol':1e-6,'maxls':50}))
    successful=[v for v in opts if v.success and np.isfinite(v.fun)]
    opt=min(successful or opts,key=lambda x:x.fun)
    attempts=len(opts)
    if not opt.success:
        attempts+=1
        opt=minimize(core,np.full(len(K),.2),jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':800,'ftol':1e-10,'maxls':50})
    if not opt.success:
        attempts+=1
        opt=minimize(core,np.maximum(opt.x,0),jac=True,method='SLSQP',bounds=bounds,options={'maxiter':1200,'ftol':1e-9})
    if not opt.success:raise RuntimeError(str(opt.message))
    f=core(opt.x,True);f.update(success=True,optimizer_attempts=attempts,iterations=int(opt.nit),upper_boundary=False)
    return f

