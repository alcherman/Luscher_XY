# Numerical systems studied with Lüscher-type admissibility conditions

Survey date: 2026-09-03.

## Executive summary

There is no single universally used meaning of “the Lüscher admissibility condition.” The numerical literature falls into several related classes:

1. **Explicit Lüscher/admissible gauge actions:** a plaquette bound is part of the integration domain, and the action is infinite or divergent at its boundary.
2. **Hard topological or continuity actions:** a nearest-neighbor angle, group distance, or plaquette angle is strictly bounded. These are direct spin-model or principal-chiral analogues of admissibility, even when the authors use the terms “constraint action,” “topological action,” or “positive-link model.”
3. **Related hard barriers:** positive-plaquette, restricted-plaquette, and bulk-preventing actions. These suppress particular lattice artifacts but generally do not satisfy the small-field bound that mathematically fixes continuum topology.
4. **Topology-fixing mechanisms that are not admissibility:** determinant ratios, monopole suppression, gradient flow, and post-processing diagnostics.

Using classes 1 and 2 as the main inclusion rule—an actual hard local restriction in the numerically sampled action or domain—I find **twelve clear simulation targets**:

| # | System | Dimension | Restriction used numerically | Representative paper |
|---:|---|---:|---|---|
| 1 | $O(2)$ quantum rotor | 1d | maximum angular increment | [INSPIRE 1355378](https://inspirehep.net/literature/1355378) |
| 2 | $O(2)$/XY model | 2d | maximum nearest-neighbor angle | [INSPIRE 1205299](https://inspirehep.net/literature/1205299) |
| 3 | $O(N)$ sigma models | 2d | maximum relative spin angle | [INSPIRE 1183821](https://inspirehep.net/literature/1183821) |
| 4 | $RP^2$ and mixed $O(3)$-$RP^2$ models | 2d | $\lvert s_x\mathbin{\cdot}s_y\rvert>C$ | [INSPIRE 397083](https://inspirehep.net/literature/397083) |
| 5 | $SU(2)$ principal chiral model | 2d | positive-link constraint | [INSPIRE 481305](https://inspirehep.net/literature/481305) |
| 6 | two-flavor massive Schwinger model | 2d | admissible $U(1)$ plaquette action | [INSPIRE 618180](https://inspirehep.net/literature/618180) |
| 7 | $U(1)$ gauge-Higgs model at $\theta=\pi$ | 2d | admissible $U(1)$ plaquette action | [INSPIRE 2807692](https://inspirehep.net/literature/2807692) |
| 8 | $\mathcal N=(2,2)$ $SU(N)$ super-Yang–Mills | 2d | Sugino admissible plaquette action | [INSPIRE 796828](https://inspirehep.net/literature/796828) |
| 9 | compact $U(1)$ pure gauge theory | 4d | hard plaquette-angle cutoff | [INSPIRE 1369107](https://inspirehep.net/literature/1369107) |
| 10 | $SU(2)$ pure Yang–Mills | 4d | hard non-Abelian plaquette action | [INSPIRE 1682468](https://inspirehep.net/literature/1682468) |
| 11 | $SU(3)$ pure gauge / quenched QCD, plus preliminary $N_f=2$ overlap work | 4d | hyperbolic/topology-conserving plaquette action | [INSPIRE 697077](https://inspirehep.net/literature/697077) |
| 12 | $SU(2)$ principal-chiral field with a Skyrme term | 3d | nearest-neighbor continuity bound preserving winding | [INSPIRE 531198](https://inspirehep.net/literature/531198) |

The number twelve is a count of system classes, not of papers. It deliberately treats the positive-link principal chiral model and spin-model constraint actions as admissibility analogues. Under a narrower convention restricted to papers explicitly calling their gauge-field bound “Lüscher admissibility,” the list would be shorter. Under a broader convention that includes any infinite local barrier, positive/restricted-plaquette $SU(2)$, $SU(4)$, and bulk-preventing $SU(5)$ actions would also enter. Those are separated below rather than silently mixed into the count.

## Foundational papers and search boundary

The historically natural **sigma-model pair** is:

- B. Berg and M. Lüscher, “Definition and Statistical Distributions of a Topological Number in the Lattice $O(3)$ Sigma Model,” *Nucl. Phys. B* **190** (1981) 412–424, [INSPIRE 164288](https://inspirehep.net/literature/164288), DOI `10.1016/0550-3213(81)90568-X`.
- M. Lüscher, “Does the Topological Susceptibility in Lattice Sigma Models Scale According to the Perturbative Renormalization Group?”, *Nucl. Phys. B* **200** (1982) 61–70, [INSPIRE 164695](https://inspirehep.net/literature/164695), DOI `10.1016/0550-3213(82)90058-X`.

Two other Lüscher papers seed partly separate gauge-theory branches:

- M. Lüscher, “Topology of Lattice Gauge Fields,” *Commun. Math. Phys.* **85** (1982) 39–48, [INSPIRE 166110](https://inspirehep.net/literature/166110), DOI `10.1007/BF02029132`.
- M. Lüscher, “Abelian chiral gauge theories on the lattice with exact gauge invariance,” *Nucl. Phys. B* **549** (1999) 295–334, [INSPIRE 479892](https://inspirehep.net/literature/479892), arXiv:hep-lat/9811032.

On the survey date, INSPIRE returned 248, 145, 359, and 308 direct citing records for those four papers, respectively, with **889 unique papers in their combined citation union**. Titles, abstracts, and metadata were screened, with full-text checks when the abstract did not say whether the restriction entered the numerical ensemble.

That census is a reproducible search boundary, not an exhaustiveness proof. Several important descendants—including constrained $RP^2$, the positive-link principal chiral model, the Skyrme simulation, and much of the Sugino literature—do not directly cite the expected historical seed. Keyword searches and citation chaining were therefore also required.

## Inclusion rule

A paper is in the main list only when the hard restriction is part of the action, update acceptance rule, or functional-integration domain used for the reported numerics. Either of the following qualifies:

- finite or zero action inside a local angle/plaquette bound and infinite action outside it;
- an action diverging at the boundary of an allowed local domain;
- an explicitly equivalent local continuity or positive-link restriction used to prevent the lattice map from crossing an ambiguity surface.

The following do **not** qualify by themselves:

- measuring the Berg–Lüscher geometric charge on configurations generated with a standard action;
- testing an admissibility inequality after gradient flow;
- fixing topology with a Wilson-kernel determinant ratio;
- suppressing monopoles without a local small-field/continuity bound;
- writing down a hard action but simulating only a smooth approximation to it.

## Main numerical literature

### 1. 1d $O(2)$ quantum rotor

- I. Bautista *et al.*, [“Measuring the Topological Susceptibility in a Fixed Sector,” INSPIRE 1355378](https://inspirehep.net/literature/1355378), arXiv:1503.06853. The paper contains Monte Carlo data for $S=0$ when every angular increment is below a fixed $\delta$ and $S=\infty$ otherwise.

W. Bietenholz *et al.*, [“Topological Lattice Actions,” INSPIRE 867576](https://inspirehep.net/literature/867576), arXiv:1009.2146, also analyze 1d $O(2)$ and $O(3)$ constraint actions, but those 1d results are analytic; the paper’s Monte Carlo component is the 2d $O(3)$ model.

### 2. 2d $O(2)$/XY model

- W. Bietenholz *et al.*, [“Topological Lattice Actions for the 2d XY Model,” INSPIRE 1205299](https://inspirehep.net/literature/1205299), arXiv:1212.0579.
- W. Bietenholz, U. Gerber and F. G. Rejón-Barrera, [“Berezinskii–Kosterlitz–Thouless transition with a constraint lattice action,” INSPIRE 1240685](https://inspirehep.net/literature/1240685), arXiv:1307.0485.
- W. Bietenholz, U. Gerber and M. Pepe, [“New insight into the Berezinskii-Kosterlitz–Thouless phase transition,” INSPIRE 1280931](https://inspirehep.net/literature/1280931), arXiv:1402.2341.
- W. Bietenholz *et al.*, [“Meron- and Semi-Vortex-Clusters as Physical Carriers of Topological Charge and Vorticity,” INSPIRE 1769183](https://inspirehep.net/literature/1769183), arXiv:1912.04112. The 2d $O(2)$ cluster-size study uses the topological/constraint-action framework.
- N. Sale, J. Giansiracusa and B. Lucini, [“Quantitative analysis of phase transitions in two-dimensional XY models using persistent homology,” INSPIRE 1926476](https://inspirehep.net/literature/1926476), arXiv:2109.10960. The numerical comparison includes the classical, topological-action, and nematic XY models.

These papers study step scaling, the helicity modulus, vortex binding and unbinding, BKT universality, cluster scaling, and persistent-homology observables.

### 3. 2d $O(N)$ sigma models

The earliest constrained numerical line includes:

- W. Bietenholz, A. Pochinsky and U.-J. Wiese, [“Meron-Cluster Simulation of the $\theta$-Vacuum in the 2-d $O(3)$-Model,” INSPIRE 395505](https://inspirehep.net/literature/395505), hep-lat/9505019. The triangular-lattice meron construction uses a maximum-angle constraint so cluster topological charges are unambiguous.
- M. Hasenbusch, [“$O(N)$ and $RP^{N-1}$ Models in Two Dimensions,” INSPIRE 397083](https://inspirehep.net/literature/397083), hep-lat/9507008. It numerically studies the hard-constraint $O(3)$ model as well as $RP^2$.
- A. Patrascioiu and E. Seiler, [“Percolation and the existence of a soft phase in the classical Heisenberg model,” arXiv:hep-th/0011199](https://arxiv.org/abs/hep-th/0011199), *J. Stat. Phys.* **106** (2002) 811. The simulated $O(3)$ action is infinite when a neighboring spin pair violates the chosen dot-product bound.
- F. Brechtefeld, [“A meron cluster solution for the sign problem of the two-dimensional $O(3)$ model,” arXiv:hep-lat/0207012](https://arxiv.org/abs/hep-lat/0207012). This continues the constrained meron-cluster simulation at $\theta=\pi$.

The modern topological/optimized-constraint line includes:

- W. Bietenholz *et al.*, [“Topological Lattice Actions,” INSPIRE 867576](https://inspirehep.net/literature/867576), arXiv:1009.2146: Monte Carlo simulations of 2d $O(3)$ using hard angular/topological actions.
- M. Bögli *et al.*, [“Non-trivial $\theta$-Vacuum Effects in the 2-d $O(3)$ Model,” INSPIRE 1080496](https://inspirehep.net/literature/1080496), arXiv:1112.1873: optimized constraint action, geometric charge, and precision step scaling at nonzero $\theta$.
- P. de Forcrand, M. Pepe and U.-J. Wiese, [“Walking near a Conformal Fixed Point: the 2-d $O(3)$ Model at $\theta$ near $\pi$ as a Test Case,” arXiv:1204.4913](https://arxiv.org/abs/1204.4913), *Phys. Rev. D* **86** (2012) 075006. The meron algorithm uses an explicit 120-degree nearest-neighbor bound.
- J. Balog *et al.*, [“Drastic Reduction of Cutoff Effects in 2-d Lattice $O(N)$ Models,” INSPIRE 1183821](https://inspirehep.net/literature/1183821), arXiv:1208.6232: numerical $O(3)$,$O(4)$,$O(8)$, with direct Monte Carlo checks for N=10,20,30,40,60, using a maximum relative-angle constraint.
- M. Bögli *et al.*, [“Study of theta-Vacua in the 2-d $O(3)$ Model,” INSPIRE 1202499](https://inspirehep.net/literature/1202499), arXiv:1211.2932: proceedings report from the same $\theta$-vacuum program.
- W. Bietenholz *et al.*, [“$O(N)$ Models with Topological Lattice Actions,” INSPIRE 1255312](https://inspirehep.net/literature/1255312), arXiv:1309.6278: proceedings overview with constraint-action data.
- I. Bautista *et al.*, [“Measuring the Topological Susceptibility in a Fixed Sector,” INSPIRE 1355378](https://inspirehep.net/literature/1355378), arXiv:1503.06853: hard-constraint 2d $O(3)$ ensembles and standard-action controls.
- W. Bietenholz, P. de Forcrand and U. Gerber, [“Topological Susceptibility from Slabs,” INSPIRE 1394397](https://inspirehep.net/literature/1394397), arXiv:1509.06433: constraint-action $O(3)$ data used to test fixed-sector/slab estimators.
- W. Bietenholz *et al.*, [“Meron- and Semi-Vortex-Clusters as Physical Carriers of Topological Charge and Vorticity,” INSPIRE 1769183](https://inspirehep.net/literature/1769183), arXiv:1912.04112: constrained 2d $O(3)$ cluster distributions and topological susceptibility.
- M. Lauk and A. Patella, [“Non-perturbative renormalization of the energy momentum tensor in the 2d $O(3)$ nonlinear sigma model,” INSPIRE 3123859](https://inspirehep.net/literature/3123859), arXiv:2602.23078, *PoS LATTICE2025* (2026) 148. After comparing three actions, all final measurements use the modified hard constraint $\cos\delta=1-1.345g_0$; the study determines EMT renormalization constants with shifted boundary conditions and a gradient-flow coupling.

Only $O(3)$ has the integer instanton charge relevant to the original Berg–Lüscher construction. The constraint-action simulations in this family cover $N=3,4,8,10,20,30,40,60$; the $N=\infty$ comparison is analytic. For $N>3$, the local angular restriction is used as an improvement or universality test rather than to protect an integer instanton charge.

### 4. 2d $RP^2$ and mixed $O(3)$-$RP^2$ models

- M. Hasenbusch, [“$O(N)$ and $RP^{N-1}$ Models in Two Dimensions,” INSPIRE 397083](https://inspirehep.net/literature/397083), hep-lat/9507008. The $RP^2$ constraint requires $\lvert s_x\mathbin{\cdot}s_y\rvert>C$, with violating links excluded; the paper reports cluster Monte Carlo and step-scaling/correlation-length data.
- C. Bonati, A. Franchi, A. Pelissetto and E. Vicari, [“Asymptotic low-temperature behavior of two-dimensional $RP^{N-1}$ models,” INSPIRE 1802502](https://inspirehep.net/literature/1802502), arXiv:2006.13061. The numerical study includes two constrained $RP^2$ models and a constrained mixed $O(3)$-$RP^2$ model.

The same 2020 paper also simulates $RP^3$, but those $N=4$ calculations use standard actions. $RP^3$ is therefore not counted as a constrained numerical target.

### 5. 2d $SU(2)$ principal chiral model: positive-link ensemble

- O. Borisenko and P. Skala, [“May vortices produce a mass gap in 2-D spin models at weak coupling?”, INSPIRE 481305](https://inspirehep.net/literature/481305), hep-lat/9812020; published as [“Vortex condensation and mass gap generation in two-dimensional principal chiral models,” INSPIRE 527971](https://inspirehep.net/literature/527971), *Phys. Rev. D* **62** (2000) 014502.

The positive-link model requires $\operatorname{Tr}(U_xU_{x+\hat\mu}^{\dagger})>0$, rejecting updates that violate the bound. The authors perform Monte Carlo simulations on lattices up to $128^2$ and compare mass gaps in the full and $Z_2$ sectors. The paper calls this the principal-chiral analogue of the positive-plaquette gauge model. It is a hard continuity constraint, although not a topology-fixing Lüscher bound in the narrow gauge-theory sense.

### 6. 2d massive two-flavor Schwinger model

- H. Fukaya and T. Onogi, [“Lattice study of the massive Schwinger model with theta term under Lüscher’s ‘admissibility’ condition,” INSPIRE 618180](https://inspirehep.net/literature/618180), hep-lat/0305004, *Phys. Rev. D* **68** (2003) 074503.
- H. Fukaya and T. Onogi, [conference version, INSPIRE 628321](https://inspirehep.net/literature/628321), hep-lat/0309104.
- H. Fukaya and T. Onogi, [“Theta vacuum effects on the chiral condensation and the eta-prime meson correlators in the two flavor massive QED(2) on the lattice,” INSPIRE 646972](https://inspirehep.net/literature/646972), hep-lat/0403024.
- H. Fukaya and T. Onogi, [“Theta vacuum effects on the pseudoscalar condensates and the eta-prime meson in 2-dimensional lattice QED,” INSPIRE 658428](https://inspirehep.net/literature/658428), hep-lat/0409021, proceedings follow-up.

These use an admissible $U(1)$ gauge action with domain-wall fermions, generate fixed sectors separately with HMC, and combine sectors to study condensates, correlators, masses, and $\theta$ dependence.

### 7. 2d $U(1)$ gauge-Higgs theory with a $\theta$ term

- S. Akiyama and Y. Kuramashi, [“Tensor renormalization group study of (1+1)-dimensional $U(1)$ gauge-Higgs model at $\theta=\pi$ with Lüscher’s admissibility condition,” INSPIRE 2807692](https://inspirehep.net/literature/2807692), arXiv:2407.10409, *JHEP* **09** (2024) 086.
- S. Akiyama and Y. Kuramashi, [proceedings follow-up, INSPIRE 2872482](https://inspirehep.net/literature/2872482), arXiv:2501.15352.

The model is evaluated with tensor renormalization rather than Markov-chain Monte Carlo. The admissible gauge action is nevertheless part of the numerical tensor formulation, so it belongs in the main list.

### 8. 2d $\mathcal N=(2,2)$ $SU(N)$ super-Yang–Mills

Sugino’s action imposes $\lVert 1-U_p\rVert<\epsilon$ and becomes singular at the boundary. Here admissibility removes lattice-artifact gauge vacua while retaining an exact lattice supercharge; topology fixing is not the primary purpose.

- H. Suzuki, [“Two-dimensional $N=(2,2)$ super Yang-Mills theory on computer,” INSPIRE 752800](https://inspirehep.net/literature/752800), arXiv:0706.1392.
- I. Kanamori, H. Suzuki and F. Sugino, [“Euclidean lattice simulation for dynamical supersymmetry breaking,” INSPIRE 767563](https://inspirehep.net/literature/767563), arXiv:0711.2099, and [longer companion, INSPIRE 767580](https://inspirehep.net/literature/767580), arXiv:0711.2132.
- I. Kanamori, [“RHMC simulation of two-dimensional $N=(2,2)$ super Yang-Mills with exact supersymmetry,” INSPIRE 795566](https://inspirehep.net/literature/795566), arXiv:0809.0655.
- I. Kanamori and H. Suzuki, [“Restoration of supersymmetry on the lattice,” INSPIRE 796828](https://inspirehep.net/literature/796828), arXiv:0809.2856.
- I. Kanamori and H. Suzuki, [“Some physics of the two-dimensional $N=(2,2)$ supersymmetric Yang-Mills theory: Lattice Monte Carlo study,” INSPIRE 802844](https://inspirehep.net/literature/802844), arXiv:0811.2851.
- I. Kanamori, [“Vacuum energy of two-dimensional $N=(2,2)$ super Yang-Mills theory,” INSPIRE 813524](https://inspirehep.net/literature/813524), arXiv:0902.2876.
- M. Hanada and I. Kanamori, [“Lattice study of two-dimensional $N=(2,2)$ super Yang-Mills at large-N,” INSPIRE 827111](https://inspirehep.net/literature/827111), arXiv:0907.4966; simulations use N=2,3,4,5.
- M. Hanada and I. Kanamori, [“Absence of sign problem in two-dimensional $N=(2,2)$ super Yang-Mills on lattice,” INSPIRE 873085](https://inspirehep.net/literature/873085), arXiv:1010.2948. Only its Sugino-action arm is an admissible ensemble; the CKKU comparison is not.

Several abstracts say only “Sugino action,” so confirming admissibility requires following the action definition into the paper or its cited formulation.

### 9. 4d compact $U(1)$ pure gauge theory

- O. Åkerlund and P. de Forcrand, [“$U(1)$ lattice gauge theory with a topological action,” INSPIRE 1369107](https://inspirehep.net/literature/1369107), arXiv:1505.02666, *JHEP* **06** (2015) 183.
- P. de Forcrand and O. Åkerlund, [“Aspects of topological actions on the lattice,” INSPIRE 1415583](https://inspirehep.net/literature/1415583), arXiv:1601.03905.

The action forbids plaquette angles beyond a cutoff. The simulations map confining, Coulomb, and monopole-free regimes and measure free energies.

### 10. 4d $SU(2)$ pure Yang–Mills

- D. Nogradi, L. Szikszai and Z. Varga, [“Non-abelian lattice gauge theory with a topological action,” INSPIRE 1682468](https://inspirehep.net/literature/1682468), arXiv:1807.05295, *JHEP* **08** (2018) 032.

This paper uses a hard non-Abelian plaquette action and tests universality against Wilson-action continuum results using the zero-temperature topological susceptibility, finite-temperature transition, and small-volume step scaling.

Earlier positive-, logarithmic-, and restricted-plaquette $SU(2)$ studies are listed below because their bounds are much looser than a topology-defining small-field condition.

### 11. 4d $SU(3)$ pure gauge / quenched QCD and preliminary dynamical overlap work

The principal numerical line studies the admissibility-derived hyperbolic or “topology-conserving” gauge action:

- S. Shcheredin *et al.*, [“Testing a topology conserving gauge action in QCD,” INSPIRE 659171](https://inspirehep.net/literature/659171), hep-lat/0409073.
- W. Bietenholz *et al.*, [“Lattice gauge actions for fixed topology,” INSPIRE 666730](https://inspirehep.net/literature/666730), hep-lat/0412017.
- K.-I. Nagai *et al.*, [“Testing topology conserving gauge actions for lattice QCD,” INSPIRE 693645](https://inspirehep.net/literature/693645), hep-lat/0509170.
- H. Fukaya *et al.*, [“Parameter dependence of the topology change and the scaling properties of the topology conserving gauge action,” INSPIRE 693741](https://inspirehep.net/literature/693741), hep-lat/0509184.
- H. Fukaya *et al.*, [“Overlap fermion with the topology conserving gauge action,” INSPIRE 695961](https://inspirehep.net/literature/695961), hep-lat/0510095. This includes preliminary two-flavor dynamical-overlap results.
- H. Fukaya *et al.*, [“Topology conserving gauge action and the overlap-Dirac operator,” INSPIRE 696316](https://inspirehep.net/literature/696316), hep-lat/0510116, *Phys. Rev. D* **73** (2006) 014503.
- W. Bietenholz *et al.*, [“Exploring topology conserving gauge actions for lattice QCD,” INSPIRE 697077](https://inspirehep.net/literature/697077), hep-lat/0511016, *JHEP* **03** (2006) 017.
- H. Fukaya, [“Lattice QCD with fixed topology,” INSPIRE 711763](https://inspirehep.net/literature/711763), hep-lat/0603008, thesis.
- F. Bruckmann *et al.*, [“Comparing topological charge definitions using topology fixing actions,” INSPIRE 820732](https://inspirehep.net/literature/820732), arXiv:0905.2849, *Eur. Phys. J. A* **43** (2010) 303.

These simulations impose a literal hard plaquette cutoff, but their practical cutoff parameters are generally much looser than the rigorous small-field bound that separates configuration space into topological sectors and prevents continuous changes of the lattice charge. They therefore report strong suppression rather than absolute elimination of topology changes. Bruckmann *et al.* also compare a determinant-ratio topology fixer; that comparison arm is not itself an admissible plaquette action.

### 12. 3d $SU(2)$ principal-chiral field with a Skyrme term

- A. J. Schramm and B. Svetitsky, [“Topology and metastability in the lattice Skyrme model,” INSPIRE 531198](https://inspirehep.net/literature/531198), hep-lat/0008003, *Phys. Rev. D* **62** (2000) 114020.
- A. J. Schramm and B. Svetitsky, [“Exact topological density in the lattice Skyrme model,” INSPIRE 533387](https://inspirehep.net/literature/533387), hep-lat/0009014, proceedings account.

This is a 3d classical lattice Skyrme model with an $SU(2)\simeq S^3$-valued field and global left-right $SU(2)$ × $SU(2)$ symmetry. The action contains the principal-chiral kinetic term and a Skyrme term. Topology preservation relies on both a divergent nearest-neighbor factor and a restricted local-update prescription; the authors choose the cutoff parameter $\alpha=0.1$ empirically. Local Metropolis Monte Carlo on a $16^3$ lattice in the $n=1$ sector is used to measure the skyrmion radius and metastability. This is an action-level fixed-winding simulation, not merely a static minimization.

## Related hard-barrier simulations: not counted in the twelve

These papers are numerically relevant if “admissibility” is being used in the broadest sense of any infinite local barrier. They are separated because the restriction does not generally supply Lüscher’s topology theorem.

### Positive, logarithmic, and restricted plaquette gauge actions

- J. Fingberg, U. M. Heller and V. Mitrjushkin, [“Scaling in the Positive Plaquette Model and Universality in $SU(2)$ Lattice Gauge Theory,” INSPIRE 375242](https://inspirehep.net/literature/375242), hep-lat/9407011: forbids negative $SU(2)$ plaquettes and studies thermodynamics, string tension, glueballs, and topological charge.
- U. M. Heller, [“$SU(2)$ lattice gauge theory with a logarithmic action,” INSPIRE 392815](https://inspirehep.net/literature/392815), hep-lat/9502009: the logarithm creates a positive-plaquette domain barrier.
- J. D. Stack and S. D. Neiman, [“Confinement by monopoles in the positive plaquette model of $SU(2)$ lattice gauge theory,” INSPIRE 415600](https://inspirehep.net/literature/415600), hep-lat/9601019.
- H. Grady, [“Monopole loop suppression and loss of confinement in restricted action $SU(2)$ lattice gauge theory,” INSPIRE 466283](https://inspirehep.net/literature/466283), hep-lat/9801016.
- A. Gavai, [“On the deconfinement transition in $SU(4)$ lattice gauge theory,” INSPIRE 584101](https://inspirehep.net/literature/584101), hep-lat/0203015. This employs a one-center-sector plaquette restriction rather than a topology-defining smoothness bound.

### Bulk-preventing gauge actions

- T. Rindlisbacher, K. Rummukainen and A. Salami, [“Bulk-preventing actions for $SU(N)$ gauge theories,” INSPIRE 2740880](https://inspirehep.net/literature/2740880), arXiv:2306.14319, *Phys. Rev. D* **108** (2023) 114511; see also the proceedings precursor [INSPIRE 1957041](https://inspirehep.net/literature/1957041), arXiv:2111.00860.

The action has an infinite barrier when a plaquette eigenvalue approaches the group cut locus. It is tested numerically in pure $SU(2)$, pure $SU(5)$, and $SU(3)$ with four Wilson-clover flavors. The authors explicitly emphasize that the action does **not** prevent topology changes, so these are not added to the main admissibility count.

## Numerical construction tests rather than equilibrium simulations

- D. Kadoh and Y. Kikukawa, [“A Numerical solution to the local cohomology problem in $U(1)$ chiral gauge theories,” INSPIRE 642914](https://inspirehep.net/literature/642914), hep-lat/0401025. The work numerically differentiates and integrates the anomaly over admissible 2d $U(1)$ gauge fields and tests locality.
- K. Nagao, [“Admissibility condition and nontrivial indices on a noncommutative torus,” INSPIRE 691389](https://inspirehep.net/literature/691389), hep-th/0509034. Random admissible matrix configurations are generated and their Ginsparg–Wilson indices computed.

These are genuine numerical uses of admissibility but not thermodynamic or equilibrium simulations of the corresponding field theory.

## Continuity-constrained soliton numerics

These calculations use the same geometric idea—exclude locally discontinuous maps so an integer degree cannot slip through the lattice—but employ numerical minimization rather than ordinary ensemble sampling.

- R. S. Ward, [“Stable topological skyrmions on the 2-D lattice,” INSPIRE 392671](https://inspirehep.net/literature/392671), hep-th/9502048: numerical minimization of 2d $O(3)$/Heisenberg baby skyrmions under an acute-neighbor-angle restriction.
- R. S. Ward, [“Hopf solitons on the lattice,” INSPIRE 699500](https://inspirehep.net/literature/699500), hep-th/0512024: 3d $S^2$ Skyrme–Faddeev/Hopf solitons with acute nearest-neighbor angles and conserved Hopf number.

The Schramm–Svetitsky Skyrme work remains in the main list because it performs Metropolis sampling of a constrained fixed-winding ensemble rather than only minimizing a classical configuration.

## Important exclusions and false positives

- **The seed papers themselves:** Berg–Lüscher define and measure a geometric topological charge, but their Monte Carlo ensemble is not generated with a hard admissibility action.
- **3d $O(4)$ in arXiv:1912.04112:** the abstract can suggest a constrained calculation, but the full text explains that the necessary bound traps the system in the ordered phase; the reported 3d $O(4)$ critical simulations use the standard action.
- **2d $RP^3$ in arXiv:2006.13061:** the constrained runs are $RP^2$ and mixed $O(3)$-$RP^2$; the $RP^3$ runs use standard actions.
- **Later Sugino/SUSY simulations:** admissibility should not be inferred merely from the name “Sugino action.” For example, [INSPIRE 1352834](https://inspirehep.net/literature/1352834), arXiv:1503.04416, explicitly avoids the admissibility restriction, and [INSPIRE 1473833](https://inspirehep.net/literature/1473833), arXiv:1607.01260, removes vacuum degeneracy without it.
- **Formulation-only work:** [INSPIRE 1088611](https://inspirehep.net/literature/1088611), arXiv:1202.2101, proposes an admissible 2d $SU(N)$ supersymmetric formulation but is not a numerical simulation paper.
- **Gradient-flow diagnostics:** Y. Tanizaki, A. Tomiya and H. Watanabe, [INSPIRE 2851198](https://inspirehep.net/literature/2851198), arXiv:2411.14812, use ordinary Wilson ensembles and test admissibility only after flow. M. Abe and O. Morikawa, [INSPIRE 2870957](https://inspirehep.net/literature/2870957), arXiv:2501.11438, likewise simulate a Wilson-type ensemble and smooth afterward.
- **Smooth approximation to a proposed hard action:** D. Banerjee *et al.*, [“An Improved Single-Plaquette Gauge Action,” INSPIRE 1409929](https://inspirehep.net/literature/1409929), arXiv:1512.04984, define a constrained action but perform the reported $SU(2)$/$SU(3)$ simulations with a finite-power smooth approximation.
- **Determinant-ratio topology fixing:** overlap-QCD simulations such as [INSPIRE 781846](https://inspirehep.net/literature/781846), arXiv:0803.3197, prohibit zero modes of the Wilson kernel with extra determinant factors. This fixes topology algorithmically but is not a Lüscher plaquette admissibility condition.
- **Monopole-free $CP^{N-1}$:** [INSPIRE 1789394](https://inspirehep.net/literature/1789394), arXiv:2003.14075, and [INSPIRE 2030100](https://inspirehep.net/literature/2030100), arXiv:2202.04614, impose zero DeGrand–Toussaint monopole charge cube by cube. This is a hard topological constraint, but it is weaker and structurally different from local small-field admissibility.
- **Homotopy-augmented fields:** M. Dall’Olio and J. M. Zapata, [INSPIRE 1866096](https://inspirehep.net/literature/1866096), arXiv:2105.14571, numerically study a 1d $U(1)$ rotor and $SO(3)$ rigid body with additional homotopy variables. The construction is explicitly an alternative to excluding non-admissible configurations.

## Unresolved historical candidate

- F. Fucito and S. Solomon, [“Instantons on the Lattice,” INSPIRE 15117](https://inspirehep.net/literature/15117), *Nucl. Phys. B* **251** (1985) 505, DOI `10.1016/0550-3213(85)90275-5`. The abstract says that modified spin- and gauge-model actions prevent tunneling under local Monte Carlo while global updates restore the relative sector weights. The accessible metadata and abstract do not identify the particular models or establish how much numerical data are reported. I therefore flag it rather than inventing additional system entries; resolving it requires the full eleven-page article.

## Completeness assessment

The report is designed to be **more complete, not formally exhaustive**. It covers:

- all titles and abstracts in the 889-record union of direct INSPIRE citations to records 164288, 164695, 166110, and 479892 as indexed on 2026-09-03;
- full-text checks for plausible numerical papers when the abstract did not specify the simulated action, with the 1985 Fucito–Solomon article explicitly left unresolved because its full text was not available in the sources inspected;
- descendant and keyword searches that recovered papers outside those citation trees;
- separate accounting for proceedings versions, numerical construction tests, classical minimization, and related but non-equivalent hard barriers.

The remaining vulnerability is bibliographic rather than conceptual: an older or cross-disciplinary paper may use an equivalent hard local restriction without the words “Lüscher,” “admissibility,” “constraint action,” or “topological action,” and without citing any of the four seed records. The positive-link principal chiral model and the $RP^2$ papers show that this happens. For that reason, the defensible conclusion is “twelve clear system classes found under the stated rule,” not “exactly twelve exist.”
