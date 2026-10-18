# Calibration harness for criteria (c) and (a) of the P1 rule

Two scripts measure how often the frozen rule gives false nuisance support, and whether it still detects planted effects, on toy cohorts (`mpaudit.sim`) and on synthetic κ laws. They were written for the design comparison of 2026-10-17 that compared fixes to the activity-selection bias of P1-GO v1.2 and led to P1-GO v1.3; they are kept as the regression check of the rule. Rerun them after any change to `mpaudit/protocol.py`, `mpaudit/decoy.py` or `mpaudit/coexpr.py`, and compare with the v1.3 reference results below. All results are toy or synthetic: they show that the code behaves as designed, not how real metaprograms behave.

- `evaluate_c.py` covers criterion (c): the protocol Δlog κ losses of `protocol.py` (null bounds, self bound, equal κ cells and depth; v1.2's activity gate).
- `evaluate_a.py` covers criterion (a): the `kappa_retained` guard of `decoy.py` (v1.3's excess over matched random sets and its entry rule; v1.2's random-set floor).

Each script runs one copy of the `mpaudit` package and writes one JSON file. That file holds `headline` (key numbers), `checks` (pass/fail against the targets below), the full summaries, per-job run times and the raw per-tumour or per-decoy results (about 2 MB for `evaluate_c.py`; keep full outputs out of the repository).

From the repository root, with the package's dependencies installed:

```
python scripts/calibration/evaluate_c.py src/mpaudit OUT_C.json [--quick] [--workers 4] [--baseline scripts/calibration/reference/c_v1.3.json]
python scripts/calibration/evaluate_a.py src/mpaudit OUT_A.json [--quick] [--workers 4] [--decoys 2]
```

`reference/` holds the v1.3 runs of this repository (`c_v1.3.json`, `a_v1.3.json`) and the v1.2 runs of commit edb2736 (`c_v1.2.json`, `a_v1.2.json`), each without the raw per-tumour rows, the module paths and the synthetic details (the headline keeps the synthetic numbers). Pass `--baseline` with a `c_*.json` of an earlier version so that the power checks (P3) are judged against it; without it they are informational.

Results come back in the order the jobs were listed, not in the order they finish, and `evaluate_a.py` sorts each group by seed and stores the seeds next to every per-seed list (`seeds`), so two runs of the same package give the same JSON. Before 2026-10-18 the per-seed lists of `evaluate_a.py` followed job completion order and carried no seed labels; `a_v1.2.json` was re-summarised from its raw rows in seed order (same values).

## How the package is loaded

- **PKG_DIR.** Give either the `mpaudit` directory itself (e.g. `src/mpaudit`, or a modified copy of it) or a directory that contains `mpaudit/`.
- **Import.** `common.load_mpaudit` drops any `mpaudit*` module that is already imported and puts PKG_DIR's parent at `sys.path[0]`, ahead of any installed copy (e.g. an editable install of `src`). It then imports `mpaudit`.
- **Import check.** `common.check_modules` asserts that every imported `mpaudit*` module lives under PKG_DIR, and the JSON records each module's path under `modules`. Workers are forked after the import, so they inherit the same modules.
- **Public API only:**

  | Script | Functions called |
  |---|---|
  | `evaluate_c.py` | `protocol.protocol_shift_kappa` (toy cohorts), `protocol.per_sample_kappa` (only to measure the realised size of injected losses), `protocol.tumour_shifts` and `protocol.protocol_shift` (synthetic levels; with `self_bound=False` when the package has that argument, i.e. from v1.3: synthetic κ levels have no cells to bootstrap, so the synthetic cases test the null bounds alone), `report.validate_rule` (frozen thresholds), `sim.simulate_cohort`, `sim.axis_magnitudes` and `data.Keys` |
  | `evaluate_a.py` | `decoy.make_decoy`, `decoy.kappa_retained`, `sim.simulate_cohort`, `data.add_qc_obs` and `report.validate_rule` |

- **Package defaults.** Package defaults are used throughout, so a package that changes a default is evaluated with its new value. This covers `n_null_sets`, `n_null_floor`, `kappa_active` (v1.2's gate; informational from v1.3) and `min_pairs`, and also `min_kappa_retained`, `min_paired_tumours` and `protocol_frac` read through `validate_rule`.
- **Synthetic set count.** The number of random sets in the synthetic cases is the default of `protocol_shift_kappa(n_null_sets=...)`.
- **`--n-null-sets N`.** Overrides the number of random sets in both the toy and the synthetic cases. It is an explicit, unregistered setting.
- **Optional hooks** (only for a package that changes the definition of criterion (c) itself):
  - `report.harness_criterion_c(per_tumour_rows_of_one_mp, rule) -> bool`
  - `report.harness_criterion_c_per_mp(per_mp, rule) -> bool array`

  Without these hooks, (c) is the frozen rule: `n_paired_tumours >= min_paired_tumours` and `frac_loss_beyond_null_kappa >= protocol_frac` (4 and 0.5).
- **`protocol.HARNESS_SYNTHETIC = False`.** A package sets this to skip the synthetic parts when its rule cannot be expressed on `(levels, null_levels)` data frames (e.g. a gate decided on held-out cells).
- **Errors.** A job that raises is recorded under `errors` and the run continues.

## Statistics used everywhere

- **Rates.** Every rate is a ratio `p = sum(k) / sum(n)`. Its Monte Carlo SE (`se`) is cluster-robust (linearised, `sqrt(K/(K-1) * sum((k_i - p n_i)^2)) / sum(n_i)`), with these clusters:

  | Where | Cluster |
  |---|---|
  | Synthetic cases | the MP, because tumours of one MP share the replicate bound |
  | Toy cases | the seed, because all sets of a seed share cells |

  Each rate is reported with `k` (events), `n` (number of tests: the denominator) and `n_clusters`.
- **Proportions over runs or MPs.** These are `P(criterion c)` and power proportions. They carry the binomial SE `sqrt(p(1-p)/n)` with `k` and `n`.
- **Definitions, per tumour–MP pair:**
  - A pair is **tested** when `loss_kappa` is not NA. NA means that the change, the null bound or (v1.3) the self bound was undefined, or (v1.2) that the activity gate failed.
  - **tested fraction** = tested pairs / all pairs.
  - **loss|tested** = pairs with `loss_kappa` True / tested pairs.
  - **gain|tested** is the same with `gain_kappa`.
  - **unconditional loss** = losses / all pairs; an untested pair counts as no loss.

## evaluate_c.py: criterion (c)

### (a) Ideal synthetic cases (`synthetic.a_per_tumour`)

Each case builds `levels` (one κ per MP and piece) and `null_levels` (one κ per random set and piece) as data frames, then calls `tumour_shifts` and `protocol_shift` with the reference protocol `fresh` and the test protocol `fixed`.

- **Pieces.** Each tumour has two pieces. In paired tumours `P*` they are fresh and fixed; in replicate tumours `R*` both are fresh.
- **Draws.** Each MP, and each of its `n_sets` random sets, draws a (piece 1, piece 2) pair from a law with correlation `rho` between the pieces. All draws are independent across tumours, MPs and sets.
- **No protocol effect** in any case.
- **Size.** 10 paired tumours. Up to 1000 MPs per case (400 with `--quick`), capped by a row budget of 12M rows (4M with `--quick`).

Laws, given as (MP, random sets):

| law | MP | random sets | meaning |
|---|---|---|---|
| `floor` | N(0.02, 0.015) | same | absent; κ mostly below v1.2's 0.05 `kappa_active` floor |
| `weak_lognormal` | lognormal(median 0.06, sdlog 0.35) | same | absent, near the floor; the scenario of `tests/test_protocol.py::test_no_gate_calibration` |
| `background` | N(0.10, 0.03) | same | absent; in v1.2 the random-set q95 set the gate's threshold |
| `high_background` | N(0.18, 0.04) | same | absent at a high background, like the calibrated toy |
| `present_weak_1.5sd` | N(0.145, 0.03) | N(0.10, 0.03) | present, weak |
| `present_3sd` | N(0.19, 0.03) | N(0.10, 0.03) | present |
| `present_strong` | N(0.30, 0.05) | N(0.10, 0.03) | clearly present |

The first four laws are "exchangeable": the MP is drawn exactly like its random sets. Each law is run at `rho` 0 and 0.5.

Paths (the names come from v1.2's activity gate; v1.3 has no gate, and every MP with replicate pairs from >= `min_pairs` tumours gets the replicate null):

| path | replicate tumours | random sets | what it tests |
|---|---|---|---|
| `random_sets` | 0 | yes | random-set null (the path most public twins will use); in v1.2 with the gate on the random-set q95 |
| `replicates_floor_gate` | 6, MP drawn as in the paired tumours | no | replicate null without random sets; in v1.2 with the floor-only gate |
| `replicates_with_sets` | 6, MP drawn as in the paired tumours | yes | replicate null with random sets available. In v1.2 an MP needed `min_pairs` tumours with a pair active in one piece and otherwise fell back to random sets; an MP exchangeable with its sets passed a q95 gate in about 5% of pieces, so it mostly fell back (see `null_types`). The rate split by the null each MP got is under `loss_given_tested_on_{replicates,random_sets}_null`. |
| `replicates_3` | 3 (the registered minimum, `min_pairs`), MP drawn as in the paired tumours | no | replicate null at its smallest df |
| `replicates_carrier_reps` | 6, MP present there (lognormal median 0.4, sdlog 0.15) | yes | the paired tumours lack the MP while the replicate tumours carry it (paired and replicate tumours not exchangeable). Informational: it violates the replicate null's stated assumption. Run only for `floor`, `weak_lognormal` and `background`. |

Reported for each law, `rho` and path: `tested_frac`, `loss_uncond`, `loss_given_tested`, `gain_given_tested`, `null_types` (MP counts per null), and `p_crit_c`.

### (b) Toy null (`toy_null`)

- **Cohort.** `simulate_cohort(n_tumours=10, n_twin=0, n_region=10, n_mal=200, n_tme_per_type=10, estimated_covariates=False, seed=s)`. Every tumour has two fresh pieces (R1, R2). There is no protocol effect at all: no stress loss, no capture bias, no selection.
- **Seeds and magnitudes.** 12 seeds (4 with `--quick`), at both the calibrated magnitudes (`cal`) and the first magnitudes (`first`, via `sim.axis_magnitudes(*FIRST_MAGNITUDES)`).
- **Gene sets per run** (11):
  - planted sets present in every tumour: cycle, hypoxia, stress;
  - planted sets present in a subset of tumours: EMT, IFN, CNA_7q;
  - 3 random 40-gene sets of unplanted generic genes;
  - the first private program of T01 and of T06.
- **Regimes:**
  - `A_random_sets`: R2 is relabelled `fixed` in all 10 tumours, so every MP gets the random-set null (no replicate tumours). Per-tumour results on this path do not depend on the other tumours.
  - `B_replicates`: R2 is relabelled `fixed` in T01–T05 only, and T06–T10 are replicate tumours (the design of the second review's null simulation). Each MP gets the replicate null if `min_pairs` tumours have a replicate pair (v1.2: an active one), otherwise random sets.
- **Calls.** One `protocol_shift_kappa(adata, keys, sets, "fresh", "fixed")` per regime, with package defaults.
- **Variants** (calibrated magnitudes, seeds 0–11, 0–1 with `--quick`; `--variant-seeds`), reported as `cal+<variant>`:
  - `test_depth_0.5`, `test_depth_0.25`: every test piece's counts (all its cells) binomially thinned to a half or a quarter of their depth, as when the test protocol is sequenced less deeply;
  - `n_mal_60`: 60 malignant cells per piece;
  - `test_cut_50`: the test pieces cut to their first 50 malignant cells.
- **Absent pairs** are tumour–set pairs where the set's program is not planted in that tumour:
  - random sets in every tumour;
  - private programs outside their owner;
  - EMT, IFN and CNA_7q in non-carrier tumours.

  **Carrier pairs** are all the others.
- **`per_tumour[mag/regime]`.** The rates above for:
  - all absent pairs, and each absent kind on its own;
  - all carrier pairs;
  - absent tested pairs active (the informational activity columns) on the reference side only (`absent_ref_side_only`) and on both sides (`absent_both_sides`), which located v1.2's bias;
  - in regime B only, absent pairs split by the null each MP got (`absent_on_replicate_null`, `absent_fallback_random_sets`).
- **`p_crit_c[mag/regime]`.** Criterion (c) per (set, seed), by category: `random`, `private_absent_everywhere` (private_T06 in regime B), `private_owner_paired`, `mixed` (EMT, IFN, CNA_7q) and `carrier_all`.
  - Regime A is evaluated on all 10 paired tumours and also on the halves T01–T05 and T06–T10, each with 5 paired tumours (`...@5_paired_tumours(halves)`).
  - `crit_c_expected_at_5pct` is the expected count if each tested tumour were a loss with probability 0.05, given the observed numbers of tested tumours.
  - `per_mp_aggregation_mismatches` counts MPs where the package's per-MP aggregation disagrees with criterion (c) recomputed from the per-tumour rows.

### (c) P(criterion c) against the number of paired tumours (`synthetic.c_p_crit_c`)

- The same synthetic generator is used with 4, 10 or 40 paired tumours.
- Laws: `weak_lognormal` (rho 0), `background` (rho 0 and 0.5) and `present_weak_1.5sd` (rho 0).
- Paths: `random_sets`, `replicates_floor_gate`, `replicates_3` and `replicates_with_sets`.
- Informational (`C2-info`): a *patchy* program, absent from a share q (0.1, 0.2, 0.3) of the pieces, each piece on its own: lognormal(median 0.4, sdlog 0.15) where present and N(0.02, 0.015) where absent, in paired and replicate tumours alike (no protocol effect), on the paths `replicates_3` and `replicates_floor_gate`, with 4 and 10 paired tumours. The replicate t bounds assume near-normal changes; these mixtures of small changes and drops to the floor are not.
- MPs per case: 2000 / 1000 / 400 for 4 / 10 / 40 paired tumours, a third of that with `--quick`.
- Each case reports:
  - `p_crit_c`: the fraction of MPs meeting criterion (c), from `protocol_shift`'s per-MP table and the frozen rule;
  - `p_powered`: the fraction with at least 4 tested tumours;
  - `p_ge3_of_20_mps_meet_c = P(Binomial(20, p_crit_c) >= 3)`, the chance that 20 independent unaffected MPs alone produce GO through (c).

### (d) Power (`power`)

- **`default`.** `simulate_cohort(seed)` for seeds 0–3 (0–1 with `--quick`), with the real protocol effect: stress is absent under the fixed protocol, plus capture bias and EMT selection. Four paired tumours (T01–T04) and two replicate tumours; all MPs use the random-set null. Reported:
  - stress tumours with a loss (`stress_tumours_detected`, out of 4 per seed);
  - stress `P(crit c)`;
  - for the other planted sets, losses among carrier and non-carrier tumours, with the list of non-carrier losses;
  - which non-stress sets met (c).
- **`inject`.** The toy-null cohort design, calibrated magnitudes, seeds 0–5 (0–1 with `--quick`), regimes A and B.
  - **The loss.** In each pseudo-fixed piece (R2 of the paired tumours), the counts of a random fraction `f` of the program's genes are permuted across the piece's malignant cells, each gene independently. This removes those genes' correlations and keeps their marginal counts.
  - **Where.** cycle in every paired tumour; EMT in the paired tumours that carry it (about half).
  - **Size.** `f` = 0.2 and 0.25 (`--inject-frac`; `--quick` runs 0.25 only). The realised κ ratio of the injected piece to the same piece before injection, measured with `per_sample_kappa`, has medians of 0.61–0.64 at f 0.2 and 0.52–0.57 at f 0.25: a κ drop of about 36–48% (EMT 36–43%).
  - **Reported** per `f{frac}/regime` and program:
    - `injected_tumours_detected`: injected tumours with `loss_kappa` True, over all injected tumours;
    - rates in carrier and non-carrier tumours (non-carriers are EMT only);
    - `p_crit_c` over seeds;
    - `kappa_ratio_injected_piece` (median, q10, q90).

### Pass criteria (`checks`)

| id | metric | target | why |
|---|---|---|---|
| C1 | max loss\|tested over every law and rho, per path (`random_sets`, `replicates_floor_gate`, `replicates_with_sets`) | <= 0.065 | 5% nominal plus Monte Carlo tolerance (SE about 0.013 at the v1.2 sizes); the review's proposed test is 0.03–0.065 |
| C1-info | the same for `replicates_carrier_reps` | <= 0.065 | informational (assumption violated) |
| C1p | tested fraction of `present_strong` | >= 0.99 | the rule must not drop clearly present MPs |
| C2 | max P(crit c) over laws x {4, 10, 40} x paths | <= 0.03 | the per-tumour 5% target gives P(>=2 of 4) = 0.014–0.03 (v1.2's floor-gate replicate path: 0.016–0.033 at 4 tumours). p_c 0.03 means P(>=3 of 20 MPs) of about 0.02. |
| C2-info | the same for the patchy laws on the replicate paths | <= 0.03 | informational: a documented gap of the replicate null (see the v1.3 results) |
| C3 | toy null loss\|tested of absent pairs, per magnitude, variant and regime | <= 0.065 | calibration on the toy |
| C4 | toy null loss\|tested of carrier pairs (C4-info for the small-piece variants `n_mal_60` and `test_cut_50`) | <= 0.02 | carriers must not be called losses; with small pieces their own κ is noisier, so that case is informational |
| C5 | toy null set x run combinations meeting (c), regime A also on halves | <= max(1, E0 + 2 sqrt(E0)), with E0 = `crit_c_expected_at_5pct` | rare-event sanity check (low power) |
| P1 | default cohort: stress tumours detected | >= 0.95 | power on the planted protocol loss |
| P2 | default cohort: stress P(crit c) | 1.0 | |
| P3 | injected loss detection per tumour, and P(crit c), per f, regime and program | >= baseline − max(2 SE, 0.05) with `--baseline`; informational without it | power must not fall below the baseline version |

A version should pass C1–C5 and P1–P2, and stay near its baseline on P3. Calibration is the priority, because false support pushes the verdict towards GO. A change that passes calibration by testing nothing shows up as a falling tested fraction (C1p) and falling P3.

## evaluate_a.py: the kappa_retained guard of criterion (a)

### Setup

- **Cohorts:**

  | cohort | simulation | sets evaluated |
  |---|---|---|
  | `default` | `simulate_cohort(seed)` | abundance, long_genes, contamination, cycle, hypoxia, EMT, IFN, CNA_7q |
  | `null` | `simulate_cohort(seed, dissociation_stress=False, tme_contamination=False, cna=False, technical_axes=False)` | cycle, hypoxia, EMT, IFN |
  | `minority` | `simulate_cohort(n_tumours=16, seed)`, 22 samples | cycle, EMT, IFN |

  `data.add_qc_obs(adata, counts_layer="counts")` adds `pct_ribo`. The default cohort has only 4–5 non-carriers, which is too few for 2 carriers among at least 8 tumours; the 16-tumour cohort gives EMT 8 non-carriers and IFN 6.
- **Arms.** The columns are sorted, as `decoy_recruitment` does:

  | arm | columns | cohorts |
  |---|---|---|
  | `technical` (registered, demo) | contamination_est, doublet_est, intronic_fraction, pct_ribo | default, null, minority |
  | `pct_ribo` | pct_ribo | default, null, minority |
  | `true_doublet` (demo `--true-doublet`) | contamination_est, doublet_score, intronic_fraction, pct_ribo | default, null |
  | `dosage` | clone_est | default |

- **Decoys.** Per (cohort, seed, arm), decoy `d` is `make_decoy(adata, keys, counts_layer="counts", strata=cols, seed=1000+d)`, for 2 decoys (`--decoys`). Seeds 0–2 (seed 0 only with `--quick`).
- **Per decoy:**
  - `kappa_retained(adata, decoy, keys, sets, counts_layer="counts", seed=1000)` with package defaults. The random-set seed 1000 is what `decoy_recruitment` uses.
  - The arm's `kappa_retained` for a seed is the mean over its decoys, as in `decoy_recruitment`.

### Metrics (`summary[cohort/arm][set]`)

- **`kappa_retained`:** mean, SE, min and max over seeds of the arm-level value, plus `per_seed`. `n_samples_whole` is the number of samples entering the median in the whole-cohort call.
- **`decomposition`** (EMT, IFN, CNA_7q). `kappa_retained` is called on each malignant sample alone, with the same full set list, so that any pool exclusions of a package match the whole-cohort call. These calls give:
  - which samples enter: carriers and non-carriers (`n_carrier_samples_entering`, `n_noncarrier_samples_entering`, `n_noncarrier_samples_total`);
  - the median over all entering samples (`kappa_retained_all_entering`) and over carriers only (`kappa_retained_carriers_only`);
  - **`noncarrier_shift_of_median`**: all-entering median minus carriers-only median, each computed per decoy and then averaged over decoys;
  - `noncarrier_ratios`: one value per entering non-carrier sample and decoy;
  - `carrier_ratio_median`.
- **`_minority[EMT|IFN]`** (minority cohort). A whole-cohort `kappa_retained` call on the restricted cohort gives `kappa_retained_as_coded`. The restricted cohort holds the first two carrier tumours by name plus every non-carrier tumour: 10 tumours for EMT, 8 for IFN, and the design for each seed is stored. Then, over every pair of carrier tumours, using the decomposition:
  - `pairs_mean_all_entering` and `pairs_max_all_entering`;
  - `pairs_mean_carriers_only`;
  - `noncarrier_shift_pairs_mean`.

### Pass criteria (`checks`; `thr` = `min_kappa_retained` = 0.5)

| id | metric | target |
|---|---|---|
| A1 | abundance, long_genes, contamination under `technical` on `default`, min over seeds | >= 0.5 |
| A2 | cycle, hypoxia, EMT, IFN under `technical` and `pct_ribo`, every cohort, max over seeds (margin = 0.5 − max is reported) | < 0.5 |
| A2-info | the same under `true_doublet` | < 0.5. Informational: cycle at about 0.50–0.67 there is the documented limit of a covariate coupled to a real program, not the floor problem. |
| A3 | CNA_7q under `dosage` on `default`, min over seeds. Dosage typing by the decoy arm also needs kappa_retained >= 0.5. | >= 0.5 |
| A4 | \|non-carrier shift of the median\| for EMT, IFN, CNA_7q under `technical` and `pct_ribo` on every cohort, max over seeds | <= 0.05 |
| A5 | minority EMT and IFN `kappa_retained` as coded, max over seeds | < 0.5 |
| A5b | minority EMT and IFN, worst pair of carrier tumours, max over seeds | < 0.5 |
| A6 | minority non-carrier shift (pairs), max over seeds | <= 0.05 |

A version should pass A1–A6. The A2 margins should not shrink, and A1/A3 must not drop below 0.5: a stricter entry rule that removes samples must not cost the technical or dosage sets their entering samples.

## Run time (4 CPUs)

| run | v1.2 (edb2736) | v1.3 (this repository) |
|---|---|---|
| `evaluate_c.py` full | 13.8 min (53 min CPU: toy null 2200 s, 24 jobs of about 92 s; synthetic 750 s; power 350 s) | 15.3 min (61 min CPU: toy null 1015 s, 24 jobs of 41–44 s, plus the variants 1582 s, 48 jobs of 25–45 s; synthetic 815 s; power 237 s). Toy null jobs are fast because each sample's random-set kappa comes from one `coexpr.set_kappas` call, the 50 bootstrap resamples of the self bound run in batches and the random-set intensities in blocks; without the variants (before 2026-10-18) the run took 8.8 min |
| `evaluate_c.py --quick` | about 4.5 min | |
| `evaluate_a.py` full | 8.0 min (27 jobs of 55–85 s) | 7.7 min (27 jobs of 55–81 s; 50 random sets per sample and set instead of 20, all kappa of a sample from one `coexpr.set_kappas` call per count matrix) |
| `evaluate_a.py --quick` | about 3.5 min | |

- Toy run time is proportional to `n_null_sets`: a version that raises it to 1000 makes the toy parts about 10 times slower. Then use `--toy-seeds`, `--parts a,c,d`, `--quick`, or run `--parts b` separately on more time.
- The synthetic row budget keeps memory and time bounded when there are more sets. It lowers the number of MPs, and the reported SEs grow accordingly.
- Run the two scripts one after the other, not at the same time: each uses `--workers` processes.

## v1.3 reference results (this repository; runs of 2026-10-18)

`reference/c_v1.3.json` and `reference/a_v1.3.json`; full runs, 4 workers, on the code of this repository. The runs of 2026-10-17 (v1.3 before equal depth, the panel on the thinned counts and the random sets drawn per tumour and MP) reproduced the design comparison's runs of the selected design value for value; those changes redraw every toy κ of (c), so the toy numbers of `c_v1.3.json` differ from that comparison's. `evaluate_a.py` is unaffected (`decoy.py`'s κ values did not change).

**evaluate_c.py: every calibration check passes** (C1–C5, C1-info and C4-info, P1, P2). The informational C2-info (a patchy program on the replicate path) fails, and P3 fails once against the v1.2 run: EMT's P(c) at f 0.25 on the random-set path, 2 of 6 runs against 5 of 6, with per-tumour detection 0.83 against 0.87 (tumours without EMT count in the denominator of (c) since v1.3; against the v1.3 run of 2026-10-17, 3 of 6, every P3 check passes).

- **Synthetic, exchangeable laws** (MP drawn like its random sets or replicate pieces), loss|tested: `random_sets` 0.020–0.054 (v1.2: 0.05–0.31), `replicates_floor_gate` 0.020–0.041, `replicates_3` 0.018–0.039, `replicates_with_sets` 0.020–0.040 (v1.2: 0.043–0.31; every MP now gets the replicate null). Present laws: 0.005–0.048. Tested fraction 1.0 everywhere (C1p). `replicates_carrier_reps`: at most 0.003 (v1.2: 0.09–0.43).
- **Synthetic P(criterion c)** for an unaffected MP: at most 0.025 on 4 paired tumours (C2; worst: `weak_lognormal`, rho 0, 3 replicate tumours; at most 0.019 on the other paths), at most 0.003 on 10, 0 on 40 (v1.2: 0.17 on 40 on the random-set path).
- **Patchy program on the replicate path (C2-info).** Absent from 10, 20 or 30% of the pieces: P(c) on 4 paired tumours 0.031, 0.045 and 0.047 with 3 replicate tumours, 0.027, 0.030 and 0.024 with 6, and at most 0.007 on 10 paired tumours. Per tumour the replicate bound stays near 5%, but all paired tumours share one bound estimated from a few heavy-tailed changes, so losses cluster within the MP. This is a documented gap of the replicate null (its t bounds assume near-normal changes), not covered by C2.
- **Toy null, pairs without the program**, loss|tested (12 seeds; every pair tested except at a quarter of the depth, where 2.5% of the pairs lose their panel):

  | case | regime A (random sets) | regime B (replicates) | carriers A / B |
  |---|---|---|---|
  | cal | 0.063 ± 0.009 (45/720) | 0.030 ± 0.012 (11/365) | 0/600, 0/295 |
  | first | 0.058 ± 0.012 (42/720) | 0.019 ± 0.007 (7/365) | 0/600, 0/295 |
  | cal, test pieces at half depth | 0.050 ± 0.011 (36/720) | 0.022 ± 0.008 (8/365) | 1/600, 0/295 |
  | cal, test pieces at a quarter of the depth | 0.048 ± 0.011 (34/702) | 0.048 ± 0.012 (17/356) | 8/600, 5/295 |
  | cal, 60 malignant cells per piece | 0.040 ± 0.009 (29/720) | 0.022 ± 0.007 (8/365) | 9/600, 1/295 |
  | cal, test pieces cut to 50 cells | 0.043 ± 0.006 (31/720) | 0.041 ± 0.008 (15/365) | 12/600, 2/295 |

  Over these 12 seeds no gene set met criterion (c) in any case or category, the 5-tumour halves of regime A included (C5). In the 24 further seeds below, one 5-tumour half did: private_T01 on T06–T10 at seed 25, where 3 of 5 tumours without the program lost (1 of 96 private-program halves; 1 of 792 set × run tests in all), within C5's tolerance. Rerun with `evaluate_c.py PKG OUT --toy-seeds 36 --parts b` (seeds 0–35; the 36-seed summary is not stored under `reference/`). The cal/A rate is higher than the 2026-10-17 run's 0.050 (36/720) because the toy κ were redrawn; 24 further seeds (12–35) gave 0.056 (81/1440), so 0.058 over 36 seeds (the design comparison's 36-seed run of the earlier code: 0.049). The depth cases are noisy by seed: at a quarter of the depth, seeds 0–5 alone gave 0.069 (A) and 0.080 (B), and seeds 6–17 0.022 and 0.017. Before equal depth (review probe of 2026-10-18, seeds 0–5), test pieces at a quarter of the depth gave 0.120 (A) and 0.27 (B), carriers 0.060 and 0.069, and an unplanted private program (A) and EMT (B) met (c); at half depth 0.071 (A).
- **Power.** Default toy (seeds 0–3): stress lost in 16 of 16 paired tumours, P(c) 4 of 4; losses of a planted set in a tumour without it: CNA_7q in T03 and EMT in T04, both seed 1 (κ fell from 0.18–0.24 to 0.05–0.06 under the real protocol), no non-stress set met (c). Injected loss (seeds 0–5): at f 0.2 (median κ ratio 0.61–0.64) cycle detected in 0.90 (A) and 0.90 (B) of the tumours, P(c) 6/6, EMT 0.47 (A, P(c) 0/6) and 0.18 (B, 0/6); at f 0.25 (κ ratio 0.52–0.57) cycle 0.98 and 0.97, EMT 0.83 (A, P(c) 2/6) and 0.53 (B, 2/6). Against v1.2 (cycle 0.92 and 0.90, 0.98 and 0.97; EMT 0.60 and 0.06, 0.87 and 0.53), per-tumour detection is within 0.13 everywhere; EMT's P(c) is lower than v1.2's (0, 0, 2, 2 of 6 against 2, 0, 5, 3) because tumours without EMT now count in the denominator of (c).

**evaluate_a.py: 58 of 60 checks pass**; the two failures are A2-info, the known coupled-covariate limit (cycle under `true_doublet`: 0.497, 0.532 and 0.530 at seeds 0, 1 and 2 on the default toy, 0.64–0.66 on the null toy). The run of 2026-10-18 gives the same values as that of 2026-10-17, every per-seed list now in seed order.

- **Technical sets under `technical`** (default toy, min over seeds): abundance 0.877, long genes 0.877, contamination 1.068 (A1).
- **Biology sets**, max over seeds: under `technical` cycle 0.031, hypoxia 0.028, EMT 0.017, IFN −0.012 on the default toy (v1.2: 0.12–0.19 on average); at most 0.036 on the minority cohort; on the null toy cycle 0.35 (seed 0; 0.07 and 0.08 at seeds 1 and 2), hypoxia 0.24, EMT 0.04, IFN 0.08. Under `pct_ribo` at most 0.001 on the default toy, 0.002 on the minority cohort and 0.16 on the null toy (A2).
- **CNA_7q under `dosage`**: 0.89–1.00 (A3).
- **Samples of tumours without the program entering**: EMT 1 of 22 and IFN 0 of 17 on the default toy, CNA_7q 0 of 14; EMT 1 of 33 and IFN 0 of 24 on the minority cohort (v1.2: 11 of 33 and 14 of 24). Shift of the median by them at most 0.012 (A4).
- **Minority case** (2 carrier tumours among 8–10): EMT −0.015 to −0.005 and IFN −0.052 to −0.021 as coded (v1.2: EMT up to 0.66); worst pair of carrier tumours 0.062 (EMT) and 0.063 (IFN); non-carrier shift at most 0.032 (A5, A5b, A6).

## v1.2 reference results (commit edb2736)

From `reference/c_v1.2.json` (full run):

- **Synthetic, random-set path, exchangeable MP**, loss|tested:

  | law | rho 0 | rho 0.5 |
  |---|---|---|
  | floor | 0.24 | 0.17 |
  | weak_lognormal | 0.31 ± 0.013 | 0.24 |
  | background | 0.18 | 0.11 |
  | high_background | 0.20 | 0.05 |

  Gain|tested is the same size, and the tested fraction is 0.04–0.12. These reproduce the review's 0.24–0.36 (lognormal) and 0.21–0.27 (normal).
- **Present MPs.** loss|tested 0.007–0.075; `present_weak` at rho 0 gives 0.075, also above the target. Tested fraction 0.64–1.0.
- **`replicates_floor_gate`:** 0.028–0.054 for exchangeable laws (0.008–0.042 for present ones). The docs' 4–5% holds for this path.
- **`replicates_with_sets`:** 0.043–0.31 for exchangeable laws. Most MPs fall back to random sets: 969–983 of 1000 for exchangeable laws at rho 0.
- **`replicates_carrier_reps`:** 0.09–0.43, which reproduces the review's 0.39–0.44 at the floor.
- **Toy null, absent pairs**, loss|tested:

  | magnitude | regime A (random sets) | regime B (replicates) |
  |---|---|---|
  | cal | 0.134 ± 0.015 (71/530; tested 0.74) | 0.096 ± 0.014 (27/280) |
  | first | 0.172 ± 0.016 (65/379; tested 0.53) | 0.103 ± 0.023 (20/195) |

  The review gave 0.136, 0.094, 0.192 and 0.104. Pairs active on the reference side only have loss|tested 0.40 (A) and 0.27–0.28 (B); pairs active on both sides 0.00–0.055.
- **Toy null, carrier pairs:** 0 losses (0/600, 0/295).
- **Toy null, criterion (c).** A random 40-gene set met (c) in 1/36 runs on 10 paired tumours at each magnitude. On 5-tumour halves it met (c) in 2/72 (cal) and 1/72 (first), and a private program outside its owner did once (1/48, first). Regime B: none. These reproduce the review's "unplanted random set meets (c)".
- **Synthetic P(crit c):**
  - random-set path at 40 paired tumours: 0.17 ± 0.019 (weak_lognormal; P(>=3 of 20) = 0.69) and 0.068 (background, rho 0);
  - 10 paired tumours: up to 0.007; 4 paired tumours: up to 0.0085;
  - floor-gate replicate path at 4 paired tumours: 0.016–0.033.
- **Power:**
  - default cohort: stress lost in 16/16 tumours, P(c) 4/4; non-carrier losses CNA_7q T03 (seeds 0 and 1) and EMT T04 (seed 1), but no non-stress set met (c);
  - injected loss at f 0.2 (κ ratio 0.62–0.64): cycle detected 0.92 (A) and 0.90 (B), P(c) 6/6; EMT 0.60 (A, P(c) 2/6) and 0.06 (B, P(c) 0/6);
  - at f 0.25 (κ ratio 0.52–0.57): cycle 0.98 and 0.97; EMT 0.87 (A, P(c) 5/6) and 0.53 (B, 3/6).
- **Checks:** C1 fails (random_sets 0.31, replicates_with_sets 0.31), C2 fails (0.17), and C3 fails for all four magnitude × regime combinations. C1 floor gate, C1p, C4, C5, P1 and P2 pass.

From `reference/a_v1.2.json` (full run, 3 seeds x 2 decoys):

- **Default cohort, technical sets under `technical`:** abundance 0.895, long_genes 0.902, contamination 1.23 (A1 passes).
- **Default cohort, biology sets:**

  | arm | cycle | hypoxia | EMT | IFN |
  |---|---|---|---|---|
  | `technical` | 0.125 | 0.147 | 0.188 | 0.145 |
  | `pct_ribo` | 0.08 | 0.11 | 0.15 | 0.11 |
  | `true_doublet` | 0.57 | | | |

- **Null cohort, cycle:** `technical` 0.17 (0.07–0.35), `true_doublet` 0.65.
- **CNA_7q:** 0.82 under `dosage` (A3 passes).
- **Non-carrier samples entering (default cohort, `technical`), per seed:**

  | set | non-carriers entering | carriers' median ratio | non-carrier ratios | shift of the median |
  |---|---|---|---|---|
  | EMT | 3–4 of 6–9 | 0.13 | 0.38–1.21 | up to +0.084 |
  | IFN | 2–3 of 5–6 | 0.12 | 0.48–0.98 | up to +0.03 |
  | CNA_7q | 1–3 of 4–5 | 0.17 | 0.21–1.37 | up to +0.03 |

  The null cohort, which has no technical axes, admits 0–3 non-carriers with low ratios.
- **Minority stress case** (2 carrier tumours among 8–10), kappa_retained as coded:

  | set | `technical` | `pct_ribo` | carriers only (pairs) |
  |---|---|---|---|
  | EMT | 0.20 / 0.13 / **0.66** (seeds 0/1/2) | up to 0.49 | 0.09–0.13 |
  | IFN | up to 0.48 | up to 0.36 | 0.11–0.14 |

  - In EMT seed 2 under `technical`, real biology passes the 0.5 guard (`tests/test_decoy_arms.py::test_minority_biology_program_stays_below_the_guard` uses this seed: v1.2 gives 0.71 there with one decoy, v1.3 0.06).
  - Worst carrier pair: 0.52–0.66.
  - Non-carrier shift: 0.23–0.39 on average over seeds, up to 0.49.
- **Checks:** A1–A3 pass. A4 fails for EMT (default and minority), and A5 (EMT `technical`), A5b and A6 fail. A2-info fails for cycle under `true_doublet` (known limit).

## Caveats

- **Synthetic cases are the ideal case.** κ per piece is an independent draw. Random sets are independent of each other and of the MP, which real sets sharing a sample's cells are not. The synthetic cases therefore test only fixes that act on `(levels, null_levels)`. Fixes that change how κ or the gate threshold are computed from cells (a held-out cell fold, excluding planted or abundance genes from the matched pool, a per-MP RNG stream) show up only in the toy parts (b) and (d) and in `evaluate_a.py`.
- **"Present" means "above its random sets"** in the synthetic laws. Between-piece heterogeneity specific to an MP (the documented gap of the random-set null) is not modelled.
- **The toy null is pseudo-paired:** both pieces are fresh and R2 is relabelled. It carries regional heterogeneity, the technical axes and background co-expression, but no fixed-protocol effects; those enter only through the `default` power cohort. Per-tumour results of regime A are independent across tumours, so the halves reuse them.
- **Toy criterion (c) counts are rare events.** C5 is a sanity check, not a calibration test; use C2 and C3 for calibration. Sets within a seed share cells, so the toy SEs are clustered by seed.
- **The injected loss is idealised.** Within-gene permutation destroys correlation but keeps means and depth relations. EMT power in regime B rests on few injected tumours (17 over 6 seeds) and is noisy. P(c) for EMT also depends on how many non-carrier tumours are tested, so a fix that tests fewer non-carriers raises EMT's P(c) without any change in per-tumour power.
- **`evaluate_a.py` does not run NMF.** Recruitment (hits, the other half of criterion (a)) is not evaluated, only the `kappa_retained` guard.
- **The decomposition is approximate.** Each per-sample call restarts the random-set seed at 1000 (v1.2: one stream over the samples; v1.3: a stream per sample and set, so this caveat applies to v1.2 only), so the samples entering can differ slightly from the whole-cohort call: compare `n_samples_whole` with the decomposition counts. Per-sample calls would also break a version that uses information across samples. The primary values (`kappa_retained`, minority `as coded`) therefore come from whole-cohort calls. The worst-pair and shift values come from the decomposition.
- **The minority "as coded" value depends on which two carriers are kept.** They are the first two by name, and carriers with two samples weigh more. The worst-pair value covers every choice.
- **Thresholds are design choices** (0.065, 0.03, 0.05), not registered settings. Report the raw values with them.
- **`--quick` is a smoke test, not a basis for decisions.** With 4 toy seeds, v1.2's cal/B_replicates absent loss|tested comes out at 0.062 (a C3 pass) against 0.096 in the full run. Decide on full runs.
- **Seeds are fixed** (the harness is deterministic for a given package), so differences between versions at the same seeds are paired comparisons.
