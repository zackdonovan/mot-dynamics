"""Paired fixed-budget search and held-out validation; run from repository root."""
import json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern,WhiteKernel,ConstantKernel
from mot.physics import *
GEOMETRY='axial6'
DIM=6
OUT=Path('results');OUT.mkdir(exist_ok=True)
def decode(x):
    if GEOMETRY=='axial6': return Controls(detuning=tuple(-.8-3.2*x[:3]),purity=tuple(.4+.6*x[3:6]))
    return Controls(detuning=tuple(-.8-3.2*x[:4]),purity=tuple([*(.4+.6*x[4:7]),2*x[7]-1]),power=(.75,)*4,geometry='oblique8')
def evaluate(x,cloud):return float(simulate(decode(x),*cloud)['score'].mean())
def save(name,data): (OUT/name).write_text(json.dumps(data,indent=2))
def main():
    start=time.time();cloud=sample_atoms(192,171); allruns=[]; chosen={}
    # Identical initial design; each algorithm receives 48 total expensive simulator calls.
    for seed in [11,23,37]:
        design=qmc.Sobol(DIM,scramble=True,seed=seed).random_base2(6)[:48]
        y0=[evaluate(x,cloud) for x in design[:16]]
        for method in ['sobol','gp_ucb']:
            xs=list(design[:16]);ys=list(y0)
            for i in range(16,48):
                if method=='sobol':x=design[i]
                else:
                    gp=GaussianProcessRegressor(kernel=ConstantKernel(1.)*Matern(length_scale=np.ones(DIM),nu=2.5)+WhiteKernel(1e-5),normalize_y=True,alpha=1e-6,optimizer=None)
                    gp.fit(np.array(xs),ys)
                    candidates=qmc.Sobol(DIM,scramble=True,seed=seed*100+i).random_base2(10)
                    mu,sd=gp.predict(candidates,return_std=True);x=candidates[np.argmax(mu+1.5*sd)]
                xs.append(x);ys.append(evaluate(x,cloud))
            best=int(np.argmax(ys));key=f'{method}_{seed}';chosen[key]=np.array(xs[best]);
            allruns.append({'method':method,'seed':seed,'x':np.array(xs).tolist(),'scores':ys,'best':best})
            print(key,max(ys),flush=True);save('search.json',allruns)
    # Model generalization separate from optimization: held-out parameter vectors.
    train=allruns[0];gp=GaussianProcessRegressor(kernel=Matern(length_scale=np.ones(DIM),nu=2.5)+WhiteKernel(1e-5),normalize_y=True,optimizer=None).fit(train['x'],train['scores'])
    tx=qmc.Sobol(DIM,scramble=True,seed=910).random_base2(5);ty=np.array([evaluate(x,cloud) for x in tx]);pred=gp.predict(tx)
    save('surrogate_test.json',{'parameter_seed':910,'n_test':32,'mae':float(abs(pred-ty).mean()),'constant_mae':float(abs(np.mean(train['scores'])-ty).mean()),'actual':ty.tolist(),'predicted':pred.tolist()})
    # No selection using evaluation samples; all six finalists reported.
    baseline=Controls() if GEOMETRY=='axial6' else Controls(detuning=(-2.,)*4,purity=(1.,1.,1.,0.),power=(.75,)*4,geometry='oblique8')
    controls={'standard':baseline,**{k:decode(x) for k,x in chosen.items()}}
    results={};testcloud=sample_atoms(4096,88291)
    for name,c in controls.items():
        result=simulate(c,*testcloud);success=result['retained'];n=len(success);p=success.mean();z=1.96;den=1+z*z/n
        center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
        results[name]={'controls':c.__dict__,'successes':int(success.sum()),'n':n,'retained_fraction':float(p),'wilson95':[center-half,center+half],'score_mean':float(result['score'].mean())}
        print(name,p,flush=True)
    save('heldout.json',results)
    # Fixed GP seed11 chosen a priori, not by test ranking; transfer perturbations.
    c=controls['gp_ucb_11'];robust={}
    for label,grad,scale in [('nominal',.12,1),('gradient_low',.096,1),('gradient_high',.144,1),('power_low',.12,.8),('power_high',.12,1.2)]:
        robust[label]={}
        for name,cc in [('standard',baseline),('gp_ucb_11',c)]:
            cfg=Controls(**{**cc.__dict__,'gradient':grad,'power':tuple(scale*p for p in cc.power)});rr=simulate(cfg,*testcloud)
            robust[label][name]=float(rr['retained'].mean())
    save('robustness.json',robust)
    small=sample_atoms(256,88292);conv={}
    for dt in [2e-4,1e-4,5e-5,2.5e-5]:
        r=simulate(c,*small,dt=dt);conv[str(dt)]={'retained':float(r['retained'].mean()),'mean_score':float(r['score'].mean()),'r':r['r'].tolist(),'v':r['v'].tolist()}
    save('convergence.json',conv);save('runtime.json',{'seconds':time.time()-start})
if __name__=='__main__':main()
