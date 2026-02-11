import os
import sys
import subprocess
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]

def run(cmd, cwd=REPO_ROOT):
    print("\n>>", " ".join(cmd))
    subprocess.check_call(cmd, cwd=str(cwd))

def exists_or_warn(path: Path, label: str):
    if not path.exists():
        print(f"[WARN] Missing {label}: {path}")
        return False
    return True

def load_yaml(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def make_paired(samples_tsv: Path, paired_out: Path):
    """
    Create a paired TSV with columns required by run_discovery.py:
    sample_id, subject_id, timepoint, response
    """
    import pandas as pd

    df = pd.read_csv(samples_tsv, sep=r"\s+", header=None, engine="python")
    if df.shape[1] < 6:
        raise ValueError(f"Expected at least 6 columns in {samples_tsv}, found {df.shape[1]}")

    # drop header row if present
    first_cell = str(df.iloc[0, 0]).strip().lower()
    if first_cell == "sample_id":
        df = df.iloc[1:].reset_index(drop=True)

    df = df.iloc[:, :6].copy()
    df.columns = ["sample_id", "subject_id", "disease", "site", "timepoint", "response"]
    out = df[["sample_id", "subject_id", "timepoint", "response"]].copy()

    paired_out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(paired_out, sep="\t", index=False)
    print(f"[OK] wrote paired file: {paired_out} ({out.shape[0]} rows)")

def main():
    import argparse
    ap = argparse.ArgumentParser(description="Run full discovery + (optional) validation pipeline.")
    ap.add_argument("--config", required=True, help="Path to YAML config (e.g., configs/full.yaml)")
    args = ap.parse_args()

    cfg = load_yaml(REPO_ROOT / args.config)

    # --- Discovery paths
    disc = cfg["discovery"]
    h5ad = REPO_ROOT / disc["h5ad"]
    samples_tsv = REPO_ROOT / disc["samples_tsv"]
    paired_out = REPO_ROOT / disc["paired_out"]
    gmt = REPO_ROOT / disc["gmt"]

    tables_dir = REPO_ROOT / cfg["outputs"]["tables_dir"]
    figs_dir = REPO_ROOT / cfg["outputs"]["figures_dir"]
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    # Preconditions
    ok = True
    ok &= exists_or_warn(h5ad, "discovery h5ad")
    ok &= exists_or_warn(samples_tsv, "discovery samples.tsv")
    ok &= exists_or_warn(gmt, "gene set gmt")
    if not ok:
        print("\n[ERROR] Missing required inputs; aborting.")
        sys.exit(2)

    # Build paired file (portable)
    make_paired(samples_tsv, paired_out)

    # Run discovery tables
    run([
        sys.executable, "scripts/run_discovery.py",
        "--h5ad", str(h5ad),
        "--paired", str(paired_out),
        "--gmt", str(gmt),
        "--outdir", str(tables_dir),
    ])

    # Figures (Discovery)
    run([sys.executable, "scripts/make_fig1_primary.py"])
    run([sys.executable, "scripts/make_fig1_means.py"])

    # --- Validation (optional)
    val = cfg.get("validation", {})
    if val.get("enabled", False):
        # Only run validation figures if required tables exist
        val_tables = REPO_ROOT / val.get("tables_dir", "results/tables/validation_gse298464")
        need1 = val_tables / "state_summary.tsv"
        need2 = val_tables / "paired_subject_deltas.tsv"

        if need1.exists():
            run([sys.executable, "scripts/make_fig1_validation_gse298464.py"])
        else:
            print(f"[WARN] Skipping figV1: missing {need1}")

        if need2.exists():
            run([sys.executable, "scripts/make_figV2_validation_gse298464_deltas.py"])
        else:
            print(f"[WARN] Skipping figV2: missing {need2}")

    print("\n[DONE] Full pipeline completed.")
    print("Tables:", tables_dir)
    print("Figures:", figs_dir)

if __name__ == "__main__":
    main()
