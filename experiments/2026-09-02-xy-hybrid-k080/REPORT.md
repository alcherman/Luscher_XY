# Hybrid update pilot at κ = 0.80

At D = 2, L = 48 and ε = π, adding one Wolff cluster to each Metropolis sweep improved the measured magnetic sampling efficiency: approximately **2.8× for |m|² and 2.0× for the Binder influence series**, relative to Metropolis alone. Energy efficiency was essentially unchanged. A fixed batch of 100 Wolff clusters remained best for magnetic observables per full runtime second, while the hybrid was more efficient for energy. These are preliminary comparisons from one cold-start chain per schedule.

![Binder estimates, effective samples per second, and runtime](comparison.png)

## Experiment and retained data

The action is $S=-\kappa\sum_{x,\mu}\cos\theta_\mu(x)$. Each method ran 10,000 fixed cycles, sequentially, with one Numba CPU thread. The schedules were:

| Name | Updates per cycle | Seed | Full generation time |
|---|---|---:|---:|
| Metropolis | One full Metropolis sweep | 2026091201 | 40.06 s |
| Hybrid | One full Metropolis sweep, then one Wolff cluster | 2026091202 | 48.02 s |
| 100 clusters | Exactly 100 Wolff clusters | 2026091203 | 236.80 s |

**Every schedule additionally makes one winding proposal at the end of each cycle.** The Metropolis proposal half-width is π. The 100-cluster count was fixed before this run using the previous observed cluster fraction of about 0.0104 at this parameter point; it was never adjusted during production. Different methods therefore have different cycle costs. The prior single-cluster L48 trajectory was not reused for these efficiency ratios.

The cold seed and all 10,000 subsequent fields, observables, random states, cluster work, and timings are saved for each arm. After inspecting the [cold-start histories](thermalization.png), the primary analysis discards the seed and cycles 1–500, retaining **cycles 501–10,000: 9,500 measurements** per arm. Energy and vortex histories settle substantially before this cut; the Metropolis magnetic relaxation is visible over roughly the first 200 cycles. See also the [full overview](overview.png) and all-observable histories in `histories/`.

## Autocorrelations and uncertainties

The Binder convention is

$$
U=1-\frac{\langle |m|^4\rangle}{2\langle |m|^2\rangle^2}.
$$

Integrated autocorrelation times use $\tau_{\mathrm{int}}=\tfrac12+\sum_{t>0}\rho(t)$, estimated on the retained series with nonperiodic FFT autocovariances and an initial-positive monotone paired sum. Values are floored at 1/2. The Binder influence series is the linearized fluctuation of this ratio of ensemble moments; its ESS diagnoses the error of the Binder estimator.

| Schedule | U ± bootstrap SE | τ energy | τ for squared magnetization | τ Binder influence | Block width / blocks |
|---|---:|---:|---:|---:|---:|
| Metropolis | 0.0445 ± 0.0348 | 6.60 | 25.05 | 12.04 | 251 / 37 |
| Hybrid | 0.0299 ± 0.0266 | 5.47 | 7.57 | 5.03 | 76 / 125 |
| 100 clusters | 0.02618 ± 0.00827 | 1.86 | 0.524 | 0.509 | 20 / 475 |

All autocorrelation times and block widths are in the schedule's own cycles. The block width is the larger of 20 and ten times the largest estimated τ among |m|², |m|⁴, Binder influence, action density and vortex density, rounded up. Binder errors use 4,000 joint contiguous-block bootstrap resamples, recomputing the ratio each time. Balanced blocks retain every measurement. The quoted widths specify the requested minimum block length; balancing distributes the remainder across the blocks. Scalar errors in `results.json` instead use equal-length batch means, omitting only the final short tail from that particular error calculation; all raw data remain saved.

The three Binder estimates are statistically compatible. The [cut and block checks](stability.png) are a useful limitation: the hybrid estimate changes to **−0.0229 ± 0.0372** when discarding 4,000 cycles. This is a shift of about two primary standard errors, using correlated subsets of the same chain. It does not establish bias, but warrants caution about the precision of the one-chain comparison. The 100-cluster estimate is less sensitive to these cuts. Increasing block lengths leaves the broad error scales similar; the largest Metropolis block choice has only nine blocks.

## Efficiency and interpretation

ESS is $n/(2\tau_{\mathrm{int}})$. The primary denominator below is the **full generation wall time**, including cold-start generation, compilation encountered during the run, measurement, HDF5 writes and checkpoint/reopen overhead.

| Schedule | Energy ESS/s | Squared magnetization ESS/s | Binder influence ESS/s |
|---|---:|---:|---:|
| Metropolis | 18.0 | 4.73 | 9.85 |
| Hybrid | 18.1 | 13.1 | 19.7 |
| 100 clusters | 10.8 | 38.3 | 39.4 |

The hybrid reduces energy autocorrelation from 6.60 to 5.47 cycles, but its extra runtime offsets that reduction. For magnetic observables, the reduction in correlation is much larger and gives a useful net gain over Metropolis.

The ranking depends on measurement cadence. Dividing by **retained generator-update time only**, including winding proposals, gives:

| Schedule | Energy ESS/update-second | Squared magnetization ESS/update-second | Binder influence ESS/update-second |
|---|---:|---:|---:|
| Metropolis | 77.2 | 20.3 | 42.3 |
| Hybrid | 75.2 | 54.3 | 81.8 |
| 100 clusters | 13.7 | 48.7 | 50.2 |

Thus the hybrid does best for magnetic observables per measured update second in this pilot, while batching 100 clusters amortizes measurement and writing costs enough to win the end-to-end comparison. This supports trying hybrid schedules, but does not identify an optimal mixture. Neither τ estimates nor efficiency ratios have uncertainty estimates here. There is one seed per method, no controlled idle-machine benchmark, and no tuning scan over proposal widths or cluster-to-sweep ratios. A short plotting process overlapped part of the cluster-only run; the wall-time figures are exploratory rather than precision hardware measurements.

## Implementation, validation, and artifacts

The runner now accepts `--update hybrid`, with `--clusters-per-step K` specifying the fixed number of clusters after each Metropolis sweep, and `--record-timing` enabling the recorded cost breakdown. Independent RNG streams and exact checkpoint replay cover both component updates and the winding proposal. Timing records do not change the sampled trajectory. The existing XY-only cluster restrictions remain in force.

The focused implementation suite passed **141 tests**:

```sh
NUMBA_NUM_THREADS=2 MPLCONFIGDIR=/private/tmp/luscher-mpl .venv/bin/python -m pytest test/test_xy.py test/test_xy_cluster.py test/test_xy_scan.py -q
```

The [saved test log](implementation-tests.txt) and [hybrid validation note](../../docs/xy-hybrid-validation.md) describe the scope. The full legacy suite was not run for this change.

The independent field audit recomputed the physical observables from **all 30,003 saved fields**, checked admissibility, canonical angles, counters, RNG key sets and timing validity, and recorded complete-file SHA256 hashes: [data-audit.json](data-audit.json). The three chains have matching nonempty source hashes in their provenance. Separate checks reproduced the retained Binder centers, bootstrap standard deviations and ESS arithmetic, and confirmed all 30 checkpoint segments completed successfully. These checks establish implementation consistency; they do not replace independent statistical replication.

Primary artifacts: [results.json](results.json), [comparison.csv](comparison.csv), [stability.json](stability.json), [autocorrelation plot](autocorrelations.png), [manifest.json](manifest.json), [status.json](status.json), [analysis provenance](analysis-provenance.json), and [PDF figure](comparison.pdf). Full raw HDF5 trajectories are in `chains/`; the scripts and unfiltered logs are preserved alongside them. Previous experiments were left untouched.

From the repository root, disk-only reproduction of the analysis is:

```sh
.venv/bin/python experiments/2026-09-02-xy-hybrid-k080/analysis.py --histories
.venv/bin/python experiments/2026-09-02-xy-hybrid-k080/analysis.py --cut 501
.venv/bin/python experiments/2026-09-02-xy-hybrid-k080/audit.py
```

`run.py` specifies the production experiment and resumes existing chains; it does not overwrite completed trajectories to create independent replicas.
