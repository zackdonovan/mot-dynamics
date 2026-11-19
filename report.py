"""Render only recorded metrics; no manual result transcription."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('results');lines=['# recorded simulation results','','These are synthetic finite-time retention outcomes, not experimental capture probabilities. Intervals describe atom sampling only. All optimizer finalists are reported; no test-set winner selection.','']
fig,ax=plt.subplots(1,2,figsize=(12,4))
for idx,(title,path) in enumerate([('six beams',root),('eight beams',root/'eight_beam')]):
    data=json.loads((path/'heldout.json').read_text());search=json.loads((path/'search.json').read_text());sur=json.loads((path/'surrogate_test.json').read_text());rob=json.loads((path/'robustness.json').read_text());con=json.loads((path/'convergence.json').read_text())
    lines+=['## '+title,'','| setting | retained / 4096 | fraction | 95% interval |','|---|---:|---:|---:|']
    for k,v in data.items():lines.append(f"| {k} | {v['successes']} | {v['retained_fraction']:.3%} | {v['wilson95'][0]:.3%}–{v['wilson95'][1]:.3%} |")
    lines+=['',f"Surrogate held-out MAE: {sur['mae']:.5f}; constant predictor MAE: {sur['constant_mae']:.5f}. Tested on 32 new parameter settings, using the training atom cloud so this isolates parameter-space interpolation rather than atom generalization.",'','| perturbation | standard | GP seed 11 |','|---|---:|---:|']
    for k,v in rob.items():lines.append(f"| {k} | {v['standard']:.3%} | {v['gp_ucb_11']:.3%} |")
    lines+=['','| step (µs) | retained fraction (256 new atoms) | smooth objective |','|---|---:|---:|']
    for k,v in con.items():lines.append(f"| {float(k)*1e6:g} | {v['retained']:.3%} | {v['mean_score']:.6f} |")
    ref=np.array(con['2.5e-05']['r']);cur=np.array(con['0.0001']['r']);err=np.linalg.norm(cur-ref,axis=1)
    lines+=['',f"100µs vs 25µs final-position discrepancy: median {np.median(err)*1000:.4f}mm, 95th percentile {np.percentile(err,95)*1000:.4f}mm. Trajectory differences and threshold crossings should not be hidden by aggregate scores.",'']
    confirmation=path/'timestep_confirm.json'
    if confirmation.exists():
        lines+=['### full holdout timestep confirmation','','Fixed comparison at 25µs, on the same 4096 held-out atoms. This refines the reported numerical outcome; it is not an extra optimization pass.','']
        for k,v in json.loads(confirmation.read_text()).items():
            lines.append(f"- {k}: 100µs {v['coarse_fraction']:.3%}, 25µs {v['fine_fraction']:.3%}; {v['changed_labels']} atom labels changed. Fine-retained atoms: 95th-percentile position discrepancy {v['position_error_p95_fine_retained_mm']:.6f}mm.")
        lines+=['','The larger all-atom position discrepancies primarily affect escaped trajectories. Aggregate retention is much more stable, but individual escaped paths should not be treated as numerically converged.','']
    for method,col in [('sobol','tab:blue'),('gp_ucb','tab:orange')]:
        curves=np.array([np.maximum.accumulate(x['scores']) for x in search if x['method']==method]);ax[idx].plot(np.arange(1,49),curves.mean(0),label=method,color=col);ax[idx].fill_between(np.arange(1,49),curves.min(0),curves.max(0),color=col,alpha=.15)
    ax[idx].set(title=title,xlabel='simulator evaluations',ylabel='best training objective');ax[idx].legend()
lines+=['Search plots show the mean best training score over three seeds; shading spans the seed minimum and maximum, not a confidence interval.','','## interpretation','','The learned model predicts simulator objective values from laser controls. Its predictions guide which expensive simulation to run next; the simulator itself is not a trained neural network. The evidence here is limited to three optimizer seeds and one declared cloud distribution. Gains over the standard setting include simple parameter tuning; gains over equal-budget Sobol search isolate the value of the surrogate more directly.','', 'The eight-beam system is not automatically superior: fixed total power is divided across more pairs, and the extra oblique pair is not optimized in geometry. This is a useful negative engineering finding, not evidence against all eight-beam MOT designs.','', 'Next scientifically meaningful work: recoil-diffusion trajectories, nonstationary internal dynamics, coherence-aware benchmark cases, calibrated source distributions and beam profiles, then measurements. Arbitrarily expanding the optimizer or claiming production use would not resolve these physical-model limits.']
Path('RESULTS.md').write_text('\n'.join(lines)+'\n');fig.tight_layout();fig.savefig(root/'search_comparison.png',dpi=160)
