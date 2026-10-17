# Four-dimensional admissible XY scout at zero coupling

This campaign measures the O(2) Binder cumulant at $\kappa=0$ on periodic
$L^4$ lattices with $L=8,12,16$, followed by selected $L=20$ checks.  Four
equally spaced points including the ends of the requested interval are used:

$$
\frac{\epsilon}{\pi}=\frac12,\frac23,\frac56,1.
$$

A second batch adds three equally spaced points inside the initial transition
bracket:

$$
\frac{\epsilon}{\pi}=\frac78,\frac{11}{12},\frac{23}{24}.
$$

A near-critical batch adds $\epsilon/\pi=0.84$ and $0.85$ before selecting
points for an $L=20$ finite-size check.

The resulting bracket selects $L=20$ at both endpoints,
$\epsilon/\pi=0.85$ and $7/8$.

A final two-point batch adds $\epsilon/\pi=0.86$ and $0.87$ on
$L=8,12,16$.  Their opposite finite-size flows tighten the scouting bracket
to $0.86<\epsilon_c/\pi<0.87$; no additional $L=20$ run is needed to classify
either point.

Each parameter point has two independent 10,000-cycle trajectories initialized
in winding sectors $(0,0,0,0)$ and $(1,0,0,0)$.  One update cycle consists of
one Metropolis sweep, four epsilon-aware Wolff clusters, and one global winding
proposal.  The runner checkpoints every 100 cycles and independently
recomputes all observables from every saved field before closure.

The Binder estimator is

$$
U_L=1-\frac{\langle|m|^4\rangle}{2\langle|m|^2\rangle^2}.
$$

Errors use a stratified contiguous-block bootstrap across the two starts, with
per-chain blocks at least ten times the largest monitored integrated
autocorrelation time.
