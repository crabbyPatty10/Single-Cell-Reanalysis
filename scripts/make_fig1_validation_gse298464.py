#!/usr/bin/env python3
"""
Validation figures for GSE298464.

Outputs:
- results/figures/figV1_gse298464_mito_ros_by_group.png
- results/figures/figV2_gse298464_mito_ros_delta.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


VAL_DIR = Path("results/tables/validation_gse298464")
OUTDIR = Path("results/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

SEED = 7  # deterministic jitter


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def _jitter(x: float, n: int, scale: float = 0.06, seed: int = SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return x + rng.normal(0.0, scale, size=n)


def _normalize_response(s: pd.Series) -> pd.Series:
    s2 = s.astype(str).str.strip()
    # common aliases
    s2 = s2.replace(
        {
            "NR": "Non_Remission",
            "R": "Remission",
            "Non-remission": "Non_Remission",
            "Non Remission": "Non_Remission",
            "Nonremission": "Non_Remission",
            "Remission ": "Remission",
        }
    )
    return s2


def _normalize_timepoint(s: pd.Series) -> pd.Series:
    s2 = s.astype(str).str.strip().str.lower()
    s2 = s2.replace(
        {
            "baseline": "Pre",
            "pre": "Pre",
            "pretreatment": "Pre",
            "before": "Pre",
            "post": "Post",
            "after": "Post",
            "week14": "Post",
            "week 14": "Post",
        }
    )
    # capitalize Pre/Post
    s2 = s2.replace({"pre": "Pre", "post": "Post"})
    return s2


def _proxy_legend(ax):
    ax.scatter([], [], marker="o", label="Non_Remission")
    ax.scatter([], [], marker="s", label="Remission")


def make_v1():
    """
    V1: %MITO_ROS-high by response and timepoint (NR Pre, NR Post, R Pre, R Post)
    Reads: results/tables/validation_gse298464/state_summary.tsv
    """
    path = VAL_DIR / "state_summary.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")

    df = pd.read_csv(path, sep="\t")

    col_response = _find_col(df, ["response"])
    col_timepoint = _find_col(df, ["timepoint"])
    col_metric = _find_col(df, ["pct_high_MITO_ROS", "pct_high_mito_ros"]) or _find_col(
        df, ["pct_high_mito_ros"]
    )

    if col_response is None or col_timepoint is None or col_metric is None:
        raise ValueError(
            "state_summary.tsv must contain columns for response, timepoint, and pct_high_MITO_ROS."
        )

    df = df.copy()
    df["response"] = _normalize_response(df[col_response])
    df["timepoint"] = _normalize_timepoint(df[col_timepoint])
    df["metric"] = pd.to_numeric(df[col_metric], errors="coerce")

    # keep only groups we expect
    df = df[df["response"].isin(["Non_Remission", "Remission"])]
    df = df[df["timepoint"].isin(["Pre", "Post"])]
    df = df.dropna(subset=["metric"])

    order = [
        ("Non_Remission", "Pre"),
        ("Non_Remission", "Post"),
        ("Remission", "Pre"),
        ("Remission", "Post"),
    ]

    data = []
    labels = []
    for r, t in order:
        vals = df[(df["response"] == r) & (df["timepoint"] == t)]["metric"].to_numpy()
        data.append(vals)
        short_r = "NR" if r == "Non_Remission" else "R"
        labels.append(f"{short_r} {t} (n={len(vals)})")

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.boxplot(
        data,
        tick_labels=labels,
        showfliers=False,  # no hollow circles
        widths=0.6,
    )

    # overlay points (jittered)
    marker_map = {"Non_Remission": "o", "Remission": "s"}
    color_map = {"Non_Remission": "tab:blue", "Remission": "tab:orange"}

    for i, (r, t) in enumerate(order, start=1):
        vals = df[(df["response"] == r) & (df["timepoint"] == t)]["metric"].to_numpy()
        if len(vals) == 0:
            continue
        xs = _jitter(float(i), len(vals), seed=SEED + i)
        ax.scatter(
            xs,
            vals,
            s=60,
            marker=marker_map[r],
            color=color_map[r],
            alpha=0.95,
            zorder=3,
        )

    _proxy_legend(ax)
    ax.legend(loc="upper right", frameon=True)

    ax.set_title("GSE298464 validation: % MITO_ROS-high by response and timepoint", fontsize=18, pad=16)
    ax.set_ylabel("% MITO_ROS-high (myeloid)", fontsize=14)

    plt.tight_layout()
    outpath = OUTDIR / "figV1_gse298464_mito_ros_by_group.png"
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Wrote: {outpath}")


def make_v2():
    """
    V2: paired delta (post - baseline) for %MITO_ROS-high, by response.
    Reads:
      - results/tables/validation_gse298464/paired_subject_deltas.tsv
      - results/tables/validation_gse298464/state_summary.tsv (for response mapping by subject_id, if needed)
    """
    dpath = VAL_DIR / "paired_subject_deltas.tsv"
    if not dpath.exists():
        raise FileNotFoundError(f"Missing {dpath}")

    d = pd.read_csv(dpath, sep="\t")

    col_subj = _find_col(d, ["subject_id"])
    col_delta = _find_col(d, ["delta_pct_high_MITO_ROS", "delta_pct_high_mito_ros"])

    if col_subj is None or col_delta is None:
        raise ValueError("paired_subject_deltas.tsv must contain subject_id and delta_pct_high_MITO_ROS.")

    d = d.copy()
    d["subject_id"] = d[col_subj].astype(str)
    d["delta"] = pd.to_numeric(d[col_delta], errors="coerce")

    # attach response if not already present
    col_resp = _find_col(d, ["response"])
    if col_resp is None:
        spath = VAL_DIR / "state_summary.tsv"
        if spath.exists():
            s = pd.read_csv(spath, sep="\t")
            s_col_subj = _find_col(s, ["subject_id"])
            s_col_resp = _find_col(s, ["response"])
            if s_col_subj is not None and s_col_resp is not None:
                s = s.copy()
                s["subject_id"] = s[s_col_subj].astype(str)
                s["response"] = _normalize_response(s[s_col_resp])
                # one response per subject (mode)
                resp_map = (
                    s.dropna(subset=["response"])
                    .groupby("subject_id")["response"]
                    .agg(lambda x: x.mode().iat[0] if len(x.mode()) else x.iloc[0])
                    .reset_index()
                )
                d = d.merge(resp_map, on="subject_id", how="left")
            else:
                d["response"] = np.nan
        else:
            d["response"] = np.nan
    else:
        d["response"] = _normalize_response(d[col_resp])

    d = d.dropna(subset=["delta"])
    d = d[d["response"].isin(["Non_Remission", "Remission"])]

    groups = ["Non_Remission", "Remission"]
    data = []
    labels = []
    for g in groups:
        vals = d[d["response"] == g]["delta"].to_numpy()
        data.append(vals)
        labels.append(f"{g.replace('_',' ')} (n={len(vals)})")

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.boxplot(
        data,
        tick_labels=labels,
        showfliers=False,  # no hollow circles
        widths=0.6,
    )

    # overlay points
    marker_map = {"Non_Remission": "o", "Remission": "s"}
    color_map = {"Non_Remission": "tab:blue", "Remission": "tab:orange"}

    for i, g in enumerate(groups, start=1):
        vals = d[d["response"] == g]["delta"].to_numpy()
        if len(vals) == 0:
            continue
        xs = _jitter(float(i), len(vals), seed=SEED + 100 + i)
        ax.scatter(
            xs,
            vals,
            s=60,
            marker=marker_map[g],
            color=color_map[g],
            alpha=0.95,
            zorder=3,
        )

    ax.axhline(0, linewidth=1)
    ax.set_title("GSE298464 validation: MITO_ROS state change by response", fontsize=18, pad=16)
    ax.text(
        0.5,
        0.96,
        "Paired delta = post − baseline",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=14,
    )
    ax.set_ylabel("Δ(Post − Pre) % MITO_ROS-high (myeloid)", fontsize=14)

    plt.tight_layout()
    outpath = OUTDIR / "figV2_gse298464_mito_ros_delta.png"
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Wrote: {outpath}")


if __name__ == "__main__":
    make_v1()
    make_v2()
