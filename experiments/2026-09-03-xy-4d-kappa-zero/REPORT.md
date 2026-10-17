# Four-dimensional zero-coupling scout

The admissible XY model was simulated at $\kappa=0$ on periodic $L^4$
lattices with $L=8,12,16$.  The original four-point grid was equally spaced on
the closed requested interval,

$$
\frac{\epsilon}{\pi}=\frac12,\frac23,\frac56,1.
$$

Three additional points divided the original transition bracket into four
equal cells:

$$
\frac{\epsilon}{\pi}=\frac78,\frac{11}{12},\frac{23}{24}.
$$

The ordered side of the transition was then resolved at
$\epsilon/\pi=0.84$ and $0.85$.  The resulting bracket motivated an $L=20$
check at both $0.85$ and $7/8$.

![Binder cumulant versus epsilon over pi](binder_vs_epsilon.png)

The final refined grid brackets the transition as

$$
0.86 < \frac{\epsilon_c}{\pi} < 0.87.
$$

At $0.86$, the Binder cumulant rises from $0.3919(32)$ at $L=8$ to
$0.48143(64)$ at $L=16$, flowing toward the ordered value $1/2$.  At $0.87$,
it falls from $0.0392(159)$ to a value consistent with zero at $L=16$.  These
opposite trends classify both new points without an additional $L=20$ run.
Linear interpolation over the sharp, unresolved drop would not be a controlled
crossing estimate.

| $\epsilon/\pi$ | $U_8$ | $U_{12}$ | $U_{16}$ | $U_{20}$ |
|---:|---:|---:|---:|---:|
| $1/2$ | 0.49998267(29) | 0.499995956(69) | 0.499998549(28) | -- |
| $2/3$ | 0.49994437(76) | 0.49998656(19) | 0.499995252(78) | -- |
| $5/6$ | 0.498808(16) | 0.4997402(36) | 0.4999160(12) | -- |
| $0.84$ | 0.498045(28) | 0.499575(6) | 0.499863(2) | -- |
| $0.85$ | 0.493687(108) | 0.498571(23) | 0.499540(8) | 0.499803(3) |
| $0.86$ | 0.3919(32) | 0.4547(17) | 0.48143(64) | -- |
| $0.87$ | 0.0392(159) | 0.0123(340) | -0.0069(266) | -- |
| $7/8$ | 0.0475(150) | 0.0346(198) | 0.0248(216) | 0.0005(366) |
| $11/12$ | -0.0089(107) | -0.0036(131) | 0.0059(125) | -- |
| $23/24$ | 0.0051(81) | 0.0183(79) | -0.0044(80) | -- |
| $1$ | 0.0044(70) | -0.0049(99) | -0.0055(77) | -- |

Each point combines two independent 10,000-cycle chains initialized in winding
sectors $(0,0,0,0)$ and $(1,0,0,0)$.  Each update cycle contains one Metropolis
sweep, four epsilon-aware Wolff clusters, and one global winding proposal.  The
first 2,001 stored rows were excluded.  Cold-start histories reach stationary
bands well before that cut; moving the cut to 1,001 or 4,001 and doubling the
block widths leaves the conclusions unchanged.

Errors are standard deviations from 4,000 stratified contiguous-block bootstrap
draws.  Each start is resampled separately, the raw moments are pooled, and the
Binder ratio is then formed.  Per-chain blocks are at least ten times the
largest integrated autocorrelation time among $|m|^2$, $|m|^4$, the Binder
influence function, link cosine, vortex density, and squared winding.  The
largest monitored $\tau_{\rm int}$ was 22.92 cycles for $|m|^4$ at
$\epsilon/\pi=7/8$ and $L=20$.  That estimate retained 34 and 46 conservative
blocks in its two starts.

All 70 trajectories passed an independent field-level recomputation of every
observable and strict admissibility checks.  Verified compaction reduced the
campaign from 257.67 GiB to 9.48 GiB while retaining every observable for all
10,001 cycles and exact fields/RNG states at 100-cycle checkpoints, analysis
cut candidates, and endpoints.
