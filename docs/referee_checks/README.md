# Referee toy simulations

Small numerical checks written by the simulated referees during the design panel
([../03_design_panel.md](../03_design_panel.md)). They are sanity checks, not validated
benchmarks, but the TESSERA proposal relies on their conclusions, so they are kept here
to be re-run.

| Script | What it checks | Result (saved output) |
|---|---|---|
| `ref1_check.py` | (1) Data thinning with a single moment-matched NB size is anti-conservative under mixtures. (2) A point-null "tumour expresses program m" test is not specific: a tumour carrying an unrelated program that shares genes with m is declared to express m. | Re-run on 2026-10-02: false membership 5% / 14% / 72% / 100% at 0 / 2 / 5 / 10 of 40 shared genes (1,000 cells, NB size 3, fold change 2.5 in 30% of cells). |
| `ref2_sim.py` | Type I error of a cross-fold co-expression test under non-NB marginals, and membership specificity under weaker effects. | `ref2_sim1.out` (type I about 0.02–0.07 with exact or moment-matched size, 0.83–0.99 with Poisson splitting); `ref2_sim2.out` (false membership about 20–50% at 4–8 of 40 shared genes). |
| `ref2_invarmp_evalue.py`, `ref2_run_one.py` | Validity and power of the split-likelihood-ratio e-value proposed in the invariance concept. | `ref2_invarmp_evalue*.out`: the e-value as specified has no power (log E far below 0 even under H1). The proposal therefore dropped it. |

Run with Python 3 and NumPy/SciPy, for example `python3 ref1_check.py`, or
`python3 ref2_run_one.py H0 0 0 1 50` (arguments: label, delta, kappa_err, seed, reps).
