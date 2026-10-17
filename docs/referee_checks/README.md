# Referee toy simulations

Small numerical checks written by the simulated referees during the design panel
([../03_design_panel.md](../03_design_panel.md)), plus one check added after the review of
2026-10-16 (`specificity_fix_check.py`). They are sanity checks, not validated benchmarks. The
TESSERA proposal ([../02_TESSERA_proposal.md](../02_TESSERA_proposal.md)) relies on their
conclusions, so every number quoted here is in a saved `.out` file in this folder. All outputs
were re-run or completed on 2026-10-16.

| Script | Question | Saved output | Result |
|---|---|---|---|
| `ref1_check.py` | (1) Is NB thinning with one moment-matched size safe under mixtures? (2) Is a point-null "tumour expresses m" test specific? | `ref1_check.out` | (1) No: fold covariance +0.125 > 0. (2) No: "membership" in 5% / 14% / 72% / 100% of tumours at 0 / 2 / 5 / 10 of 40 genes shared with an unrelated program. kappa(m) > 0 is true for k >= 2, so this is a specificity failure of the estimand, not a type I error. |
| `ref2_sim.py sim1` | Type I error of the cross-fold co-expression test when marginals are not NB | `ref2_sim1.out` | 0.027–0.067 with exact or moment-matched NB size; 0.83–1.00 with Poisson splitting. |
| `ref2_sim.py sim2` | The same specificity failure with a weaker program, by number of cells; does an empirical null help? | `ref2_sim2.out`, `ref2_sim2_seed1.out`, `ref2_sim2_seed3.out` | At 8 of 40 shared genes: 28% / 50% / 89% at 300 / 1,000 / 3,000 cells (3 seeds pooled). The empirical null leaves it unchanged. |
| `ref2_invarmp_evalue.py`, `ref2_run_one.py` | Validity and power of InvarMP's thinned split-likelihood-ratio e-value | `ref2_invarmp_evalue_runs.out`, `ref2_invarmp_evalue.out` | The e-value as specified never fired under a null but had no power at a 1.65-fold shift (it detects 2.7-fold). The plug-in variant has power but fires under the null when gene-wise capture varies. |
| `specificity_fix_check.py` | Does the proposal's fix remove the ref1 failure? The fix is kappa conditional on the other programs, a matched-random null band and a floor of 0.05. | `specificity_fix_check.out`, `specificity_fix_check_depth.out` | Yes when the overlapping program is in the dictionary: <= 1% false membership up to 30 of 40 shared genes, with power 1.00. Without conditioning, only the floor protects, and it fails once the overlap pushes kappa above 0.05 (from 20 of 40 shared genes with the true depth). |

## ref1_check.py → `ref1_check.out`

This check runs in about 30 s (seed 1, 100 replicates per row), and the re-run reproduces the
numbers quoted since 2026-10-02 exactly.

1. **Thinning with one moment-matched size.** The data are a 50/50 mixture of NB(mean 100,
   size 100) and NB(mean 1, size 0.5). The moment-matched size is 1.020, larger than the 0.5
   of the second component. Thinning with it gives that component a fold covariance of
   +0.1251 by Monte Carlo (200,000 cells), against +0.1286 in theory. Positive fold covariance
   means the thinned test is anti-conservative.
2. **Specificity of a point-null membership test.**
   - Setup: 1,000 cells, 300 genes, all with mean 2, NB size 3. The tumour carries program w
     (genes 0–39, fold change 2.5 in 30% of cells) and never expresses m. m has 40 genes, k of
     which are w genes.
   - Test: exact NB thinning into two halves and a one-sided cross-fold z-test of
     kappa(m) > 0.
   - Result: P(reject) = 0.05 / 0.14 / 0.72 / 1.00 at k = 0 / 2 / 5 / 10.

**How to read the 72% and 100%.** At k = 0, m's genes are independent of each other and of w,
so 0.05 is the type I error, and it is nominal. For k >= 2, kappa(m) > 0 is **true**: the k
shared genes really are co-expressed through w. In raw latent expression,
kappa(m) = k(k−1)/1560 × 0.355, which is 0.0005, 0.005 and 0.020 at k = 2, 5 and 10. The test
is right about its own null hypothesis. What fails is the estimand, because "m's genes co-vary
in this tumour" does not mean "this tumour expresses m". The 72% and 100% are therefore a
specificity failure of the point-null estimand, not a type I error. They should not be set
against a type I error target. They motivate changing the estimand, and the false-membership
target in 02 §6.3 applies to the new estimand.

## ref2_sim.py → `ref2_sim1.out`, `ref2_sim2.out`, `ref2_sim2_seed1.out`, `ref2_sim2_seed3.out`

Usage: `python ref2_sim.py {sim1|sim2} SEED`. The seed defaults to 1.

**sim1: type I error under non-NB marginals** (`ref2_sim1.out`, seed 1, 4.6 min).
- Setup:
  - 300 cells, 400 genes, and no cross-gene structure (H0 true).
  - Three marginals: exact NB with size 3; independent per-gene on/off states (8% of cells
    "on", 8-fold) with NB size 10; lognormal(sd 1)-Poisson.
  - Two selection engines on fold A: the top 30 loadings of PC1, or a seed gene plus the 29
    genes most correlated with it.
  - Dirichlet-multinomial thinning 50/20/30.
  - 300 replicates per row, 100 for Poisson splitting.
- The re-run reproduces the 13 previously saved rows exactly and adds the 5 missing lognormal
  rows.
- Results:
  - Exact or moment-matched size: type I 0.027–0.067.
  - Half the moment-matched size: 0.023–0.033. Mean Z is −0.28 to −0.31, so the test is
    conservative at a cost in power.
  - Poisson splitting: 0.83–1.00.
- The largest departure with moment matching is not a rejection rate. With lognormal-Poisson
  data and the PC engine, type I is 0.037 but mean Z is +0.20, which is 3.6 standard errors
  (SE 0.96/√300).

**sim2: membership specificity with a weaker program.**
- Setup:
  - The tumour expresses program P (40 genes, fold change 2.0 in 30% of cells). NB size 3 is
    known exactly. Depth is lognormal(0, 0.4) and gene means are lognormal(0, 1).
  - m consists of k genes of P plus 40 − k other genes.
  - "Empirical null": the cross-fold z is standardised by the z of 200 random 40-gene sets
    drawn from the non-P genes.
  - 100 replicates per row (50 at 3,000 cells).
- `ref2_sim2.out` is seed 2. That run reproduces the 11 previously saved rows exactly and adds
  the missing row: 3,000 cells at 8/40, 0.90 → 0.90.
- Seeds 1 and 3 are independent repeats. Pooled results over the three seeds (point null →
  after empirical null):

| shared genes | 300 cells | 1,000 cells | 3,000 cells |
|---|---|---|---|
| 0/40 | 0.06 → 0.06 | 0.08 → 0.07 | 0.05 → 0.03 |
| 2/40 | 0.04 → 0.05 | 0.06 → 0.07 | 0.03 → 0.05 |
| 4/40 | 0.10 → 0.10 | 0.11 → 0.12 | 0.17 → 0.19 |
| 8/40 | 0.28 → 0.29 | 0.50 → 0.54 | 0.89 → 0.89 |

- Monte Carlo SE of the pooled rates is 0.011–0.03.
- Individual seeds vary more than the earlier text implied. For example, 4/40 at 1,000 cells
  is 0.05, 0.18 and 0.10 in seeds 1, 2 and 3.

**What sim2 shows.**
- As in ref1, kappa(m) > 0 is true for k >= 2, so the rejections measure the specificity of
  the estimand.
- False membership grows with the number of cells (0.28 → 0.50 → 0.89 at 8/40), so
  point-null membership inflates with depth.
- Standardising by random gene sets does not help (8/40: 0.27 → 0.29 and 0.47 → 0.52 in
  `ref2_sim2.out`). Random non-P sets do not carry the overlap, so they cannot absorb it.

## ref2_invarmp_evalue.py, ref2_run_one.py → `ref2_invarmp_evalue_runs.out`, `ref2_invarmp_evalue.out`

**Setup.**
- One tumour with 3,000 cells in environment 0 and 800 in environment a; 300 genes; 3
  programs; NB size 3.
- Exact NB thinning X = X1 + X2 with the true size, then KL-NMF.
- The numerator is the X1 unrestricted fit evaluated on X2 env-a cells. The denominator is the
  restricted fit, in one of two variants:
  - "spec", as specified in the concept: per-cell usage refitted on X2;
  - "plug": X1 usage on both sides.
- Under H1, program 1 has a deviation of exp(δ) on 20 of its 40 genes in environment a.
- "kappa SD" is a gene-wise capture error of that size on the log scale. It is shared by all
  programs, and there is no program-specific shift.
- Each row reports log E under a Poisson working likelihood and under the NB likelihood.

**`ref2_invarmp_evalue_runs.out`** was produced by `ref2_run_one.py LABEL DELTA KAPPA_ERR SEED
REPS` from the design panel. Its seeds and replicate counts were not recorded. 03 says 6
replicates per row, and seed 1 with 6 replicates does not reproduce the file. Mean log E
(Poisson / NB), with P(E >= 20) in brackets:

| Row | spec | plug |
|---|---|---|
| H0, δ = 0 | −8,859 / −2,652 (0.00) | −184 / −94 (0.00) |
| H1, 1.65-fold (δ = 0.5) | −5,618 / −1,776 (0.00) | +3,690 / +910 (1.00) |
| H1, 2.72-fold (δ = 1.0) | +7,973 / +1,521 (1.00) | +20,828 / +4,354 (1.00) |
| H0, kappa SD 0.2 | −5,548 / −1,543 (0.00) | +3,987 / +1,417 (1.00) |

**`ref2_invarmp_evalue.out`** is the script's own loop (seed 7, 20 replicates per row,
11 min). The earlier file held only its first row. The re-run reproduces that row exactly and
completes the other four:

| Row | spec | plug |
|---|---|---|
| H0, δ = 0 | −8,716 / −2,654 (0.00) | −133 / −77 (0.00) |
| H1, 1.65-fold | −5,808 / −1,737 (0.00) | +3,841 / +906 (1.00) |
| H1, 2.72-fold | +9,178 / +1,589 (1.00) | +21,005 / +4,357 (1.00) |
| H0, kappa SD 0.1 | −7,961 / −2,413 (0.00) | +887 / +333 (1.00) |
| H0, kappa SD 0.2 | −5,563 / −1,488 (0.00) | +3,944 / +1,484 (1.00) |

**Conclusion.** The earlier statement that "the e-value as specified has no power" was wrong.
What the outputs show:
- **The specified e-value** never reached E >= 20 under any null, including misestimated
  capture, so it is valid in these runs. It detects a 2.7-fold shift on 20 genes (log E
  +7,973 Poisson / +1,521 NB) but not a 1.65-fold shift. Under H0, log E sits between −8,700
  and −8,900 (Poisson) and near −2,650 (NB), because the denominator refits per-cell usage on X2 and the
  numerator cannot. That incidental-parameter penalty makes the e-value very conservative and
  underpowered at moderate effects.
- **The plug-in e-value** restores power at 1.65-fold. However, it fires with P(E >= 20) = 1.00
  under the null whenever gene-wise capture is misestimated, already at SD 0.1. It is invalid
  under nuisance misspecification.

Dropping universal-inference e-values (02, risk 15) is still right. The reason is that the
powered variant is invalid under misspecification and the valid variant is underpowered, not
that the e-value has no power.

## specificity_fix_check.py → `specificity_fix_check.out`, `specificity_fix_check_depth.out`

This check is new (2026-10-16). It tests the fix proposed in 02 §1.2 definition 1 and §2.4:
kappa_tm is computed after removing the activities of every other program, and tumour t
expresses m if H0: kappa_tm <= max(0.05, q95 of kappa over matched random gene sets in t) is
rejected (one-sided, cell-bootstrap SE).

**Data.**
- As in ref1: 1,000 cells, 300 genes, NB size 3, program w = genes 0–39, and m with 40 genes,
  k of them shared with w. k = 0, 2, 5, 10, 20, plus 30 for non-member tumours.
- An active program multiplies its genes' means by 2.5 in an independent 30% of cells. When w
  and m are both active on a shared gene, the effects add.
- Scenarios:
  - "w only": t does not express m, so every rejection is a false membership;
  - "m only" and "m and w": t expresses m, so rejections measure power.
- 200 replicates per row. The main run takes 17 min with 3 workers, and the depth run 6.3 min
  with 2.
- Two baseline settings:
  - "flat": exactly ref1, every gene with mean 2;
  - "varied": gene means 2 × lognormal(sd 0.5).

**Arms.**
- **A, naive.** ref1's point-null cross-fold z-test.
- **B2, w unknown.** No conditioning. kappa is the mean off-diagonal of
  `mpaudit.coexpr.latent_corr` (CS-CORE moment estimator, size factor = total UMI). The null
  band is the q95 of 200 matched sets from `mpaudit.coexpr.matched_random_sets`, the floor is
  0.05, and 50 bootstrap draws give the SE.
- **B1, w known.** The dictionary {w, m} has the true loadings and the true baseline (a
  perfect dictionary). Each cell's activities come from NNLS on its expression proportions.
  w's activity enters the CS-CORE mean model as a covariate, so kappa is a latent partial
  correlation given w. The null band uses the same conditioning.

**Results, main run** (`specificity_fix_check.out`, size factor = total UMI). False membership
of non-member tumours ("w only"):

| k shared (of 40) | 0 | 2 | 5 | 10 | 20 | 30 |
|---|---|---|---|---|---|---|
| kappa(m), raw latent expression, theory | 0 | 0.0005 | 0.005 | 0.020 | 0.086 | 0.198 |
| A, naive (flat / varied) | 0.04 / 0.06 | 0.07 / 0.10 | 0.74 / 0.75 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 |
| B2, null band + floor, no conditioning | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 1.00 / 1.00 |
| B1, conditional on w + null band + floor | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.01 |
| mean kappa estimate, B2 (flat) | +0.013 | +0.004 | −0.005 | −0.009 | +0.027 | +0.121 |
| mean kappa estimate, B1 (flat) | −0.010 | −0.009 | −0.009 | −0.008 | −0.009 | −0.010 |

- **Null band.** The mean null q95 is −0.007 to 0.018 in non-member tumours, so the floor of
  0.05 binds in 100% of those replicates. The null band never decides anything there.
- **Depth normalisation in this toy.** An active w raises a cell's total UMI by about 20%
  (it is 40 of 300 genes). Normalising by total UMI pulls shared and non-shared genes of m
  apart, which is why B2's kappa (0.027 at k = 20) is well below the raw value (0.086). In real
  data with thousands of genes, a 40-gene program usually moves total UMI much less.
- **Power** at fold change 2.5 (kappa about 0.27, far above the floor) is 1.00 for every arm
  in "m only" and "m and w", with one exception. B2 in the flat setting with w present has
  power 0.52 / 0.36 / 0.05 / 0.27 / 1.00 at k = 0 / 2 / 5 / 10 / 20. With flat baselines, the
  genes of m and w are the only genes with raised means, so mean-matched random sets for m's
  genes are drawn from w's genes. The null q95 is then 0.23–0.28, as large as kappa(m) itself.
  This is a toy artefact of a 300-gene universe, but it shows that mean matching can draw
  controls from another active program. B1 conditions the null sets on w as well and keeps
  power at 1.00.
- **Attenuation from conditioning.** Conditioning on w's estimated activity costs m some
  signal. In "m only" at k = 0, mean kappa is 0.273 for B2 and 0.242 for B1. Cells where m is
  active have lower proportions of w's genes, so NNLS gives them lower w activity. The
  attenuation shrinks as k grows and is harmless at this effect size.

**Results, true depth** (`specificity_fix_check_depth.out`, "w only", size factor = the true
depth, which is the same for every cell). This removes the compositional shrinkage. B2's
kappa then matches theory (0.087 at k = 20, 0.199 at k = 30), and B2 false membership becomes:

| k shared | 0 | 2 | 5 | 10 | 20 | 30 |
|---|---|---|---|---|---|---|
| B2 (flat / varied) | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | **0.98 / 1.00** | 1.00 / 1.00 |
| B1 (flat / varied) | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.01 |

B1's mean kappa is between −0.001 and +0.018, so the −0.01 offset in the main run came from
total-UMI normalisation.

**What it shows.**
- **Conditioning on the overlapping program removes the failure.** When w is in the
  dictionary, B1 has false membership <= 1% (0 of 200 at every k <= 20, and 0–2 of 200 at
  k = 30) with power 1.00. It removes the overlap signal itself: kappa is about 0 at every k.
- **The within-tumour null band contributes nothing here.** Its q95 stays below the floor in
  non-member tumours.
- **Without conditioning, the floor kappa_min = 0.05 is the only protection.** It holds only
  while the overlap-induced kappa stays below 0.05. With a program of this strength, theory
  puts the crossing at k >= 16 of 40. The simulation fails at 20 of 40 with the true depth,
  and at 30 of 40 with total-UMI normalisation.

**Caveats.**
- The dictionary is perfect.
- Arm B2 assumes w is modelled nowhere. In 02 it could also be removed as a private program.
- There is a single tumour, so no BH step.
- The membership fold is all 1,000 cells. 02's 40/30/30 split would leave 300.
- The gene universe has only 300 genes, so random sets contain about 5 w genes.
- Power was checked only at an effect far above the floor. The power cost of the floor for
  weak programs was not simulated.

## Implications for 02

- **02 §1.3, the specificity bullet.** It says conditioning plus an effect-size null "removes
  that failure". That holds only for the conditioning part, and only when the overlapping
  program is in the dictionary (B1, <= 1%). The empirical null of `ref2_sim2.out` and the
  matched-random null band in the new check do not help. The floor helps only while the
  overlap kappa is below 0.05.
- **02 §6.3, the specificity row.** The ref1 figures (72% at 5/40, 100% at 10/40) are rejections
  of a null that is false. Comparing them with a <= 5% false-membership target contrasts two
  estimands; it is not a type I error comparison. The comparable figures for the proposed
  estimand are the B1 and B2 rows above.

## Running

Run from this folder, with Python 3, NumPy and SciPy. `ref2_run_one.py` loads
`ref2_invarmp_evalue.py` by relative path, and `specificity_fix_check.py` imports `mpaudit`
(`pip install -e .` at the repository root).

```
python ref1_check.py > ref1_check.out                                  # 30 s
python ref2_sim.py sim1 1 > ref2_sim1.out                              # 4.6 min
python ref2_sim.py sim2 2 > ref2_sim2.out                              # 8.8 min; seeds 1 and 3 -> ref2_sim2_seed{1,3}.out
python ref2_invarmp_evalue.py > ref2_invarmp_evalue.out                # 11 min (seed 7)
python ref2_run_one.py H0_kappaSD0.2 0 0.2 SEED REPS                   # one row of ref2_invarmp_evalue_runs.out
python specificity_fix_check.py --workers 3 > specificity_fix_check.out   # 17 min
python specificity_fix_check.py --size-factor depth --scenarios "w only" --workers 2 > specificity_fix_check_depth.out   # 6 min
```
