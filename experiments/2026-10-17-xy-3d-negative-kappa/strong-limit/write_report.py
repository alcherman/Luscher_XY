from pathlib import Path
import json
P=Path(__file__).resolve().parent
a=json.loads((P/'asymptote-assessment.json').read_text())
lines=['3D negative-kappa asymptote reconnaissance, 2026-10-17','',
'Headline: The attained finite-box responses are consistent with a narrow competing-response region approaching epsilon/pi=1/2 on the scale 1/abs(kappa). They do not establish an equilibrium thermodynamic transition curve, its uniqueness, or a numerical asymptote estimate. The scan discovers vector suppression AND recovery, coexisting with strong staggered-nematic response. Initial-condition trapping remains material.','',
'Sampling and scope:',
'96 N6 trajectories at kappa=-4,-8,-16; 60 matching selected N8 trajectories plus eight fixed-q=.60 N6/N8 controls at kappa=-16; 48 N8 trajectories at kappa=-32,-64 including eight fixed-q=.60 controls. All retain fields, native/higher-harmonic observables, RNG states, exact checkpoints and full logs.',
'The first 164 trajectories use 8000 cycles, with 2000 warmup. The final 48 use 16000 cycles, with 4000 warmup. Four separately seeded initial-state classes: uniform zero winding, unit-x winding, canted checkerboard, vortex-rich parity-locked Z4. The last is not fully disordered: its initial staggered nematic and fourth-harmonic powers equal 1.',
'Each cycle uses one single-thread local SiteUpdate sweep, four signed reflection-cluster updates and one winding-update attempt. These are one-slot runs: no replica exchanges occur. Numerical native/model kernels are unchanged. Full independent closure audits establish data/update bookkeeping and source consistency, not mixing.',
'At most 8 simulation workers, preserving half the 18 cores when combined with one analysis process. The previously paused exponent campaign remains paused.','',
'Observables and statistics:',
'Maximum first-harmonic Fourier power is evaluated per configuration then averaged, sampling every 50 cycles; this detects nonzero-momentum spin response missed by uniform magnetization. It carries a finite-volume extreme-value noise floor and is not a critical estimator. This differs from the maximum of the averaged spectrum. Uniform/staggered magnetization, staggered second-harmonic power, fourth-harmonic power, vorticity, mean link cosine and winding histories are retained.',
'All initial-state classes remain identifiable. Conditional block bootstrap uses 1000 draws and conservative blocks based on measured postcut/late-half autocorrelation times. Low/insufficient-block flags and split-half diagnostics are retained. No precise transition intervals are extracted. Plotted ranges span initial classes, not confidence intervals.',
'Fixed-u grid coordinates approach 1/2 by construction. Only agreement of measured responses at fixed u and the independent fixed-q controls can supply evidence. No response threshold is claimed universal.','',
'Attained counts:',f"{a['trajectories']} trajectories; {a['low_block_count']} low-block flags; {a['insufficient_blocks']} insufficient-block flags.",f"Generation+closure-audit total {a['generation_plus_audit_worker_seconds']/3600:.4f} worker-hours, not elapsed wall time or CPU utilization. Actual full-file bytes {a['actual_bytes']} ({a['actual_bytes']/2**30:.2f}GiB).",'',
'Selected N8 responses (ranges across uniform/canted/Z4 starts; winding starts separately in figures):',
'kappa  u  epsilon/pi  max-spin-power range  staggered-nematic range  vortex-density range']
for d in a['N8_scaled_comparison']:
 lines.append(f"{d['kappa']:4d} {d['u']:4.2f} {d['q']:.8f} {d['three_start_vector_range']} {d['three_start_nematic_range']} {d['three_start_vortex_range']}")
lines+=['','Energy-sector sanity check:',
'An explicit admissible bipartite product box bounds the partition function from below. For eta=1/(12abs(kappa)) and c=cos(epsilon-4eta)+(log(pi/eta)+1)/(3abs(kappa)), equilibrium P(mean-link-cosine>=c)<=exp(-N^3). See energy_sector_check.py and asymptote-assessment.json for the derivation, parameters and actual frame counts. Some winding-start histories, the kappa=-64 uniform start at u=1, and the kappa=-64 fixed-q=.60 Z4 start persist in these exponentially atypical sectors; the latter two each have 12000/12000 postcut frames above threshold. The three non-unit classes are not a selected equilibrated reference. This is overwhelming evidence against interpreting those histories as stationary Gibbs samples. No independence of consecutive frames is assumed: for a stationary history, the expected fraction in the event is bounded by the same probability. Failure to hit the event does not certify equilibrium or order of the complementary trajectories.',
'', 'Heuristic interpretation, not a theorem about phase order:',
'A plaquette vortex requires the sum of four principal bond angles to reach 2pi. The strict cutoff therefore forbids vortices at epsilon<=pi/2. Near epsilon=pi/2+delta, a saturated-plaquette estimate gives an extra action of approximately 4abs(kappa)*delta relative to the bond-energy infimum. This motivates testing u=abs(kappa)*(epsilon-pi/2). At fixed positive delta, large negative kappa can suppress vortices again, so the observed spin-response recovery is physically plausible. Neither a local core-cost argument nor finite-box agreement proves infinite-volume order, excludes other ordered phases, or identifies a unique continuation of the original transition.',
'', 'Interpretation limits:',
'Near-unit staggered nematic response means the spin-suppressed region must not be labeled fully disordered. The remaining winding/other-start differences and few effective blocks prevent a defensible equilibrium asymptote fit. A next scientific step would require improved equilibration and larger-volume scaling, rather than quoting more digits from these small boxes. No additional production batch is silently queued.',
'', 'Artifacts:',
'negative-kappa-asymptote-response-map.png/pdf/svg: physical coordinates, separate initial states.',
'strong-negative-scaled-response.png/pdf/svg: response alignment and winding-start failures.',
'N8-response-table.csv and asymptote-assessment.json: numerical data, sources, energy-event counts and qualifications.',
'All predecessor data, failed logs and reviewed snapshots remain preserved.']
(P/'negative-kappa-asymptote-assessment-2026-10-17.txt').write_text('\n'.join(lines)+'\n')
