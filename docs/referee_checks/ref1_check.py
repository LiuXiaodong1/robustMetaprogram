import numpy as np
rng=np.random.default_rng(1)
# (1) Corollary 1 counterexample: component-specific NB size
w=np.array([.5,.5]); mu=np.array([100.,1.]); b=np.array([100.,0.5])
m=(w*mu).sum(); EV=(w*(mu**2/b)).sum(); Vmu=(w*mu**2).sum()-m**2
bmarg=m**2/(EV+Vmu); print("b_marg(moments)=%.3f vs component sizes"%bmarg, b)
# Monte Carlo: thin with b~=bmarg, check within-component-2 covariance of folds
def thin(X,bt,eps=(0.5,0.5)):
    out=np.zeros((len(X),2))
    for i,x in enumerate(X):
        if x==0: continue
        p=rng.dirichlet(np.array(eps)*bt); out[i]=rng.multinomial(x,p)
    return out
n=200000
X2=rng.negative_binomial(b[1], b[1]/(b[1]+mu[1]), size=n)
F=thin(X2,bmarg)
emp=np.cov(F[:,0],F[:,1])[0,1]
eA=eB=.5; th=eA*eB*mu[1]**2*(bmarg-b[1])/(b[1]*(bmarg+1))
print("component-2 fold cov: emp=%.4f theory=%.4f (positive => anti-conservative)"%(emp,th))
# (2) kappa>0 triviality: tumour with program w (no m), m shares k genes with w
G=300; n=1000; nb=3
FC=2.5
def sim(k_shared, reps=200):
    rej=0
    for r in range(reps):
        wgenes=np.arange(0,40); mgenes=np.concatenate([np.arange(40-k_shared,40), np.arange(100,100+40-k_shared)])
        act=rng.random(n)<0.3
        lam=np.full((n,G),2.0); lam[np.ix_(act,wgenes)]*=FC
        X=rng.negative_binomial(nb, nb/(nb+lam))
        # exact thinning A/C with true b (0.5/0.5 for simplicity)
        A=np.zeros_like(X); 
        for g in range(G):
            pa=rng.beta(0.5*nb,0.5*nb,size=n); A[:,g]=rng.binomial(X[:,g],pa)
        C=X-A
        def res(Y,eps):
            mu_=Y.mean(0,keepdims=True); return (Y-mu_)/np.sqrt(mu_+mu_**2/(eps*nb))
        RA=res(A,.5); RC=res(C,.5)
        mm=np.zeros(G); mm[mgenes]=1; mm/=np.linalg.norm(mm)
        h=(RA@mm)*(RC@mm)-(RA*RC)@(mm**2)
        Z=np.sqrt(n)*h.mean()/h.std()
        rej+= Z>1.645
    return rej/reps
for k in [0,2,5,10]:
    print("m shares %2d/40 genes with an unrelated program w: P(reject 'tumour expresses m') = %.2f"%(k,sim(k,100)))
