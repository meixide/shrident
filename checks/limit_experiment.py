# Limit experiment under exact identification: Z ~ N(0, Omega), T = max(Q(Z),0),
# sqrt(n)(psi_EB - psi*) -> sum_j w_j(T) Z_j.  Compare with Gaussian claim and sandwich coverage.
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'apps'))
from inference import convex_weights
rng=np.random.default_rng(1)
def run(J,sig,rho,R=200000):
    Om=rho*np.outer(sig,sig)+(1-rho)*np.diag(sig**2)
    Z=rng.multivariate_normal(np.zeros(J),Om,size=R)
    s2=sig**2; iu=np.triu_indices(J,1)
    Q=(((Z[:,:,None]-Z[:,None,:])**2-(s2[:,None]+s2[None,:]))[:,iu[0],iu[1]]).mean(1)/2
    T=np.maximum(Q,0)
    W=1/(s2[None,:]+T[:,None]); W/=W.sum(1,keepdims=True)
    L=(W*Z).sum(1)
    Vs=np.einsum('rj,jk,rk->r',W,Om,W)          # sandwich at tau2-hat (conditional on T)
    cov=np.mean(np.abs(L)<=1.645*np.sqrt(Vs))
    w0=1/s2; w0/=w0.sum(); v0=w0@Om@w0           # Gaussian claim with tau0^2=0 weights
    k=np.mean(((L-L.mean())/L.std())**4)
    gls=1/np.sum(np.linalg.inv(Om)); best=s2.min()
    wc=convex_weights(Om); cvx=wc@Om@wc              # minimum of w'Omega w over the simplex
    return dict(P_T0=np.mean(T==0).round(3),var_limit=L.var().round(3),var_claimed=round(v0,3),
                kurtosis=round(k,3),sandwich_cov90=round(cov,4),var_GLS=round(gls,3),var_best_single=round(best,3),var_convex=round(float(cvx),3),var_precision_tau0=round(v0,3))
for J,rho in [(5,0.0),(5,0.5),(5,0.9),(10,0.5)]:
    print(J,rho,run(J,np.linspace(0.6,2.0,J),rho))
print('two-estimator example v=(1,4),rho=.9:')
s=np.array([1.,2.]);Om=.9*np.outer(s,s)+.1*np.diag(s**2);w=(1/s**2)/(1/s**2).sum()
print(' IVW var',round(w@Om@w,3),' best',1.0,' GLS',round(1/np.sum(np.linalg.inv(Om)),3))
