"""Command line: ``mpaudit p1`` runs the pilot audit; ``mpaudit simulate`` writes a toy cohort."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _load_config(path: str | None) -> dict:
    if not path:
        return {}
    text = Path(path).read_text()
    if path.endswith((".yaml", ".yml")):
        import yaml
        return yaml.safe_load(text) or {}
    return json.loads(text)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="mpaudit", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("p1", help="run the P1 audit on an .h5ad cohort")
    p1.add_argument("h5ad")
    p1.add_argument("--config", help="YAML or JSON config (see configs/p1_example.yaml)")
    p1.add_argument("--out", default="results/p1")
    p1.add_argument("--n-jobs", type=int)

    sim = sub.add_parser("simulate", help="write the toy cohort (and its cytoband) for a demo run")
    sim.add_argument("--out", default="results/sim")
    sim.add_argument("--seed", type=int, default=0)
    sim.add_argument("--n-tumours", type=int, default=10)

    args = ap.parse_args(argv)
    import anndata as ad

    if args.cmd == "p1":
        from .pipeline import run_p1
        cfg = _load_config(args.config)
        if args.n_jobs:
            cfg["n_jobs"] = args.n_jobs
        run_p1(ad.read_h5ad(args.h5ad), cfg, args.out)
    elif args.cmd == "simulate":
        from .sim import simulate_cohort
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        adata, truth = simulate_cohort(n_tumours=args.n_tumours, seed=args.seed)
        adata.write_h5ad(out / "sim_cohort.h5ad")
        truth["cytoband"].to_csv(out / "cytoband.txt", sep="\t", header=False, index=False)
        from .data import write_gmt
        write_gmt(truth["sets"], out / "truth.gmt")
        (out / "truth_kinds.json").write_text(json.dumps(truth["kinds"], indent=2))
        print(f"wrote {out}/sim_cohort.h5ad ({adata.n_obs} cells x {adata.n_vars} genes), cytoband.txt, truth.gmt")


if __name__ == "__main__":
    main()
