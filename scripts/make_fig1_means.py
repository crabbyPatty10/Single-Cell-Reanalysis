import os
import pandas as pd
import matplotlib.pyplot as plt

def boxplot_by_group(deltas: pd.DataFrame, colname: str, outname: str, ylabel: str, title: str) -> None:
    nr = deltas.loc[deltas["resp_group"] == "Non_Remission", colname].dropna()
    rr = deltas.loc[deltas["resp_group"] == "Remission", colname].dropna()
    print(f"{colname}: Non_Remission n={len(nr)} median={float(nr.median()):.4f} | Remission n={len(rr)} median={float(rr.median()):.4f}")

    plt.figure()
    plt.boxplot([nr, rr], tick_labels=["Non_Remission", "Remission"], showfliers=True)
    plt.axhline(0, linewidth=1)
    plt.ylabel(ylabel)
    plt.title(title)
    outpath = os.path.join("results", "figures", outname)
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print("Saved:", outpath)

def main() -> None:
    # Ensure we run relative to repo root even if executed elsewhere
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(repo_root)
    print("CWD:", os.getcwd())

    deltas_path = os.path.join("results", "tables", "paired_subject_deltas.tsv")
    deltas = pd.read_csv(deltas_path, sep="\t")
    os.makedirs(os.path.join("results", "figures"), exist_ok=True)

    boxplot_by_group(
        deltas,
        "delta_mean_ROS",
        "fig1c_delta_mean_ROS.png",
        "Delta(Post - Pre) mean ROS module score",
        "TAURUS discovery: Mean ROS change by response",
    )

    boxplot_by_group(
        deltas,
        "delta_mean_OXPHOS",
        "fig1d_delta_mean_OXPHOS.png",
        "Delta(Post - Pre) mean OXPHOS module score",
        "TAURUS discovery: Mean OXPHOS change by response",
    )

if __name__ == "__main__":
    main()
