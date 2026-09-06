"""Run a prespecified QC/graph grid into separate result directories."""

import argparse
import copy
import subprocess
import sys
from pathlib import Path

import yaml


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--mt-thresholds", type=float, nargs="+", default=[15, 20, 25])
    p.add_argument("--spatial-k", type=int, nargs="+", default=[4, 6, 8])
    p.add_argument("--execute", action="store_true", help="Omit to generate reviewable configs only")
    args = p.parse_args()
    original = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    plan_dir = Path("work/sensitivity")
    plan_dir.mkdir(parents=True, exist_ok=True)
    for mt in args.mt_thresholds:
        for k in args.spatial_k:
            if not 0 < mt <= 100 or k < 1:
                raise ValueError("Invalid QC threshold or neighborhood size")
            cfg = copy.deepcopy(original)
            cfg["qc"]["max_mt_pct"] = mt
            cfg["analysis"]["spatial_k"] = k
            cfg["output_dir"] = original["output_dir"] + f"_mt{mt:g}_k{k}"
            path = plan_dir / f"mt{mt:g}_k{k}.yaml"
            path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
            print(path)
            if args.execute:
                subprocess.run([sys.executable, "-m", "bcflow.cli", "run", "--config", str(path)],
                               check=True)


if __name__ == "__main__":
    main()

