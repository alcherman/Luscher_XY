# Independent review of the XY runtime forecast

An independent agent reviewed the forecast, its implementation, the companion
JSON and calibration job records. The arithmetic was reproduced: 36 audited
trajectories, 7.272 GiB of data, 1.3642 generation worker-hours and 0.1090 audit
worker-hours; 324 baseline trajectories including calibration, 54 selected
size checks and 342 new trajectories. The proposed new d=3,4 work projects to
18.653 worker-hours and 104.563 GiB under the stated scaling assumptions.

The review accepted cautious first-pass feasibility before Friday. It did
not support a claim that complete validated +/-0.02 curves are guaranteed.
All timing pilots use epsilon/pi=1; constrained-slice cost remains unmeasured.
The 12–36-hour figure is a planning allowance, not a measured full-campaign
runtime interval. Unmeasured d=2 additions and adaptive estimator work prevent
a measured upper bound.

The report was revised to state prominently that pilot timings are from the
ordinary XY endpoint, take no credit for overnight 18-worker operation, and
describe the 22 d=2 slices as baseline candidates pending endpoint validation.
The original internal-storage scenario left only 14.690 GiB above the reserve;
that capacity limitation was removed by the subsequently verified SSD route.
SSD writing and audit throughput still need actual measurement. This review
does not validate the transition estimator or its infinite-volume uncertainty.
