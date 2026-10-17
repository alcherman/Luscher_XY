"""Source-backed summary of the attained negative-kappa response scan."""
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
a=json.loads((P/'analysis.json').read_text());b=json.loads((P/'response-bands.json').read_text())
rows=[r for r in a['rows'] if r['point'].get('role')=='production']
lines=['3D constrained XY: measured negative-kappa continuation reconnaissance','',
'Purpose: visualize how uniform-order response continues from kappa=0 toward epsilon/pi=1/2. The scan does not approach epsilon=pi. No asymptote is imposed.', '',
'MODEL AND SAMPLING',
'Action S=-kappa*sum_positive_links cos(delta phi), with strict |principal delta phi|<epsilon. Periodic N=8 boxes in D=3.',
'At each fixed epsilon, configurations exchange among a ladder of negative couplings. Each slot cycle has one frozen Metropolis sweep, four signed-coupling reflection clusters, and the native global twist attempt. Adjacent exchange parity alternates. The signed bond rule and full product exchange invariant measure have independent exact tests; these correctness checks do not prove mixing.',
'Each of seven epsilon slices has two fresh independent ladders: zero winding and unit-x winding initial fields. All thermodynamic slots of an ensemble are correlated; they are not independent replicas for statistical inference.',
'16,000 cycles per ladder; first4,000 excluded in this descriptive analysis. All field configurations, measurements, RNG states, provenance and unfiltered logs are preserved. Joint checkpoints commit all slots at one epoch. Exact checkpoint replay and field/swap audits are tested.',
f"Main scan: {len(rows)} independent ladders, {sum(len(r['point']['kappas'])*r['point']['steps'] for r in rows):,} slot cycles, all14 closure audits passed. Analysis uses shared time blocks across all slots, and independent ensembles for bootstrap replication.", '',
'ATTAINED RESPONSE BANDS (NOT CONFIDENCE INTERVALS)',
b['definition'],
'epsilon/pi    kappa_low    kappa_high',]
for r in b['bands']:
 if r.get('role','production')!='production':continue
 lines.append(f"{r['epsilon_over_pi']:.2f}          {r.get('kappa_low','unresolved')}          {r.get('kappa_high','unresolved')}")
lines += ['', 'INTERPRETATION AND LIMITS',
'The measured finite-box response changes move toward smaller epsilon as kappa decreases. This is compatible with the proposed continuation toward epsilon/pi=1/2; it does not establish that limit, a thermodynamic critical curve, or a unique phase below it.',
'The q=.52 scan has no common stable response band. Winding trapping persists; it must not be treated as evidence of no transition or of a bound below kappa=-2.',
'All main scan ensembles trigger the conservative ladder diagnostic: minimum adjacent rate below0.05 OR fewer than2 full round trips. The q=.75 ensembles actually made9/16 round trips, q=.78 made92/81; q<=.70 made0. Do not conflate the diagnostic flag with zero trips everywhere. Some also have insufficient blocks, start disagreement or half shifts. Thus narrow interpolated Binder roots are retained as diagnostics but excluded from the headline figure. Stable responses at separated anchors do not rule out shared metastability; these are empirical bands, not rigorous equilibrium brackets.',
'Rare ordered excursions and very negative Binder values are compatible with metastability/coexistence, but do not demonstrate a first-order transition. Small uniform magnetization is not assumed to imply a massive disordered phase: nonzero Fourier modes are recorded separately.',
'No finite-size extrapolation of the negative branch has been accepted. N12 pilot failures are retained; the earlier positive/zero-kappa attained estimates are reused unchanged with their original conditional MC and descriptive sensitivity qualifications.',
'A separate two-ladder N8 q=.65 overlap validation uses .0025 central spacing, .01 outer spacing, two fresh starts,16,000 cycles and a shorter [-.40,-.25] coupling ladder. Read its completed audit and mixing diagnostics separately; no success is assumed from the design.', '',
'REPRODUCIBILITY',
'analysis.json contains all mean observables, half-series checks, tau/block diagnostics, Fourier warnings, walker round trips, bootstrap roots and root-multiplicity flags. response-bands.json records all four-half anchors and exact source IDs. Raw files are in the dedicated SSD chains symlink. Producer, audit, runner and immutable source snapshots are retained with batch manifests. No historical paused exponent jobs were resumed.',
'Plot: measured-negative-kappa-continuation.png/pdf/svg; no fitted or imposed negative asymptote.',
f"analysis SHA256: {hashlib.sha256((P/'analysis.json').read_bytes()).hexdigest()}",
f"response bands SHA256: {hashlib.sha256((P/'response-bands.json').read_bytes()).hexdigest()}"]
(P/'negative-kappa-response-scan-report.txt').write_text('\n'.join(lines)+'\n')
fine=[r for r in a['rows'] if r['point'].get('role')=='fine_transport_validation']
if fine:
 with (P/'negative-kappa-response-scan-report.txt').open('a') as f:
  f.write('\nCOMPLETED CENTRAL-OVERLAP CHECK\nBoth fresh validation ensembles passed closure audits. Minimum adjacent exchange rates: '+str([min(r['exchange_rates']) for r in fine])+'. Full walker round trips: '+str([r['walker_round_trips']['round_trips'] for r in fine])+'. Descriptive large/small uniform-component changes along walker configurations: '+str([sum(r['walker_uniform_response_changes']) for r in fine])+'. These are correlated events, not independent equilibrium phase-conversion samples.\nThe dense quarter-Binder interpolation is near kappa=-.337. Its narrow conditional bootstrap interval is NOT accepted: full transport fails, independent-start offsets remain, and one ensemble has only3 conservative blocks. Keep the broad observed band[-.40,-.30] in the headline figure. No critical-coupling precision improvement is claimed. See overlap-validation-summary.json and independent-start-response-checks.png/pdf/svg.\n')
