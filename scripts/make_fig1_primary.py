import os
import pandas as pd
import matplotlib.pyplot as plt


def main():
    # Always run relative to repo root
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(repo_root)
    os.makedirs(r"results\figures", exist_ok=True)

    deltas = pd.read_csv(r"results\tables\paired_subject_deltas.tsv", sep="\t")

    # --- Figure 1: delta MITO_ROS-high by response ---
    nr = deltas.loc[deltas["resp_group"] == "Non_Remission", "delta_pct_high_MITO_ROS"].dropna()
    rr = deltas.loc[deltas["resp_group"] == "Remission", "delta_pct_high_MITO_ROS"].dropna()

    print("Fig1 medians:",
          "Non_Remission", float(nr.median()), "|",
          "Remission", float(rr.median()))

    plt.figure()
    plt.boxplot([nr, rr], tick_labels=["Non_Remission", "Remission"], showfliers=True)
    plt.axhline(0, linewidth=1)
    plt.ylabel("Delta(Post - Pre) % MITO_ROS-high (myeloid)")
    plt.title("TAURUS discovery: MITO_ROS state change by response")
    out1 = r"results\figures\fig1_mito_ros_delta.png"
    plt.savefig(out1, dpi=200, bbox_inches="tight")
    plt.close()
    print("Saved:", out1)

    # --- Figure 1b: delta MITO_ROS vs delta fraction inflamed ---
    x = deltas["delta_frac_inflamed"]
    y = deltas["delta_pct_high_MITO_ROS"]

    plt.figure()
    plt.scatter(x, y)
    plt.axhline(0, linewidth=1)
    plt.axvline(0, linewidth=1)
    plt.xlabel("Delta(Post - Pre) fraction inflamed biopsies")
    plt.ylabel("Delta(Post - Pre) % MITO_ROS-high (myeloid)")
    plt.title("TAURUS discovery: MITO_ROS change vs inflammation change")
    out2 = r"results\figures\fig1b_mito_ros_vs_inflammation.png"
    plt.savefig(out2, dpi=200, bbox_inches="tight")
    plt.close()
    print("Saved:", out2)


if __name__ == "__main__":
    main()
