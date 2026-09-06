"""Command-line entry points; failures are explicit and never become fake results."""

import argparse
import logging
import os
from pathlib import Path

from . import __version__
from .core import config, mkdir, write_json


def main():
    parser = argparse.ArgumentParser(description="Patient-aware breast cancer multiomics")
    parser.add_argument("command", choices=["demo-data", "singlecell", "pseudobulk", "bulk",
                                            "spatial", "integrate", "figures", "run"])
    parser.add_argument("--config", default="config/demo.yaml")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()
    cfg = config(args.config)
    out = mkdir(cfg["output_dir"])
    os.environ.setdefault("MPLCONFIGDIR", str(mkdir("work/matplotlib").resolve()))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(
                            out / "pipeline.log", encoding="utf-8")])
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    from . import demo, differential, figures, integration, singlecell, spatial
    stages = {"demo-data": demo.generate, "singlecell": singlecell.run,
              "pseudobulk": differential.pseudobulk, "bulk": differential.bulk,
              "spatial": spatial.run, "integrate": integration.run, "figures": figures.run}
    selected = list(stages) if args.command == "run" else [args.command]
    if cfg["mode"] == "real" and "demo-data" in selected:
        selected.remove("demo-data")
    write_json(out / f"{args.command}_status.json", {"status": "running", "mode": cfg["mode"],
               "command": args.command, "stages": selected})
    try:
        for stage in selected:
            logging.info("Starting %s [%s]", stage, cfg["mode"])
            result = stages[stage](cfg)
            logging.info("Completed %s: %s", stage, result or "OK")
        write_json(out / f"{args.command}_status.json", {"status": "completed", "mode": cfg["mode"],
                   "command": args.command, "config_path": str(Path(args.config)), "stages": selected})
    except Exception as exc:
        logging.exception("Pipeline failed; no downstream completion is claimed")
        write_json(out / f"{args.command}_status.json", {"status": "failed", "error": str(exc),
                   "mode": cfg["mode"], "command": args.command})
        raise


if __name__ == "__main__":
    main()
