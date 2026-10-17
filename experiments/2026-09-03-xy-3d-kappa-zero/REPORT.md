# Three-dimensional zero-coupling scout

The admissible XY model was simulated at $\kappa=0$ on periodic $L^3$
lattices with $L=8,12,16$.  The epsilon grid was
$2.40,2.45,2.50,2.55,2.60$ radians; the duplicate $2.50$ in the request was
interpreted as $2.55$.

![Binder cumulant versus epsilon](binder_vs_epsilon.png)

The Binder curves reverse their volume ordering between $\epsilon=2.50$ and
$2.55$, bracketing the finite-volume crossing in that interval.  Straight-line
interpolation over this coarse cell gives crossings of 2.50799 for $L=8,12$,
2.50666 for $L=12,16$, and 2.50742 for $L=8,16$.  Their agreement is striking,
but these are scouting estimates: the 0.05 epsilon spacing and finite-volume
corrections dominate the precision of a critical-cutoff estimate.

| $\epsilon$ | $U_8$ | $U_{12}$ | $U_{16}$ |
|---:|---:|---:|---:|
| 2.40 | 0.494127(86) | 0.497454(37) | 0.498638(21) |
| 2.45 | 0.486378(225) | 0.493780(95) | 0.496768(52) |
| 2.50 | 0.439441(1276) | 0.460349(913) | 0.473385(707) |
| 2.55 | 0.242992(4784) | 0.132997(8940) | 0.048147(13912) |
| 2.60 | 0.084474(9664) | 0.022736(14952) | 0.013293(21671) |

Each point combines two independent 10,000-cycle chains initialized in winding
sectors $(0,0,0)$ and $(1,0,0)$.  Each cycle contains one Metropolis sweep,
four epsilon-aware Wolff cluster updates, and one global winding proposal.  The
first 2,001 stored rows were excluded.  The complete and early-time histories
show stationary bands before that cut, and results remain stable when the cut
is moved to 1,001 or 4,001.

Errors are standard deviations from 4,000 stratified contiguous-block bootstrap
draws.  Each independent start is resampled separately, raw moments are pooled,
and the Binder ratio is formed after pooling.  Per-chain block widths are at
least ten times the largest integrated autocorrelation time among $|m|^2$,
$|m|^4$, the Binder influence function, link cosine, vortex density, and squared
winding.  The largest monitored $\tau_{\rm int}$ was 7.56 cycles; all estimates
had at least 105 blocks per start, and doubling the block widths left the errors
stable.

All 30 chains passed an independent field-level recomputation of observables and
strict admissibility checks.  The original dense records occupied 6.1 GiB.
Verified compaction reduced them to 249 MiB while retaining every observable at
all 10,001 cycles and exact fields/RNG states at 100-cycle checkpoints, analysis
cut candidates, and endpoints.
