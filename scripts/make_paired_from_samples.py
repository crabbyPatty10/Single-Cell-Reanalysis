import pandas as pd
from pathlib import Path

IN_PATH = Path(r"data\discovery\taurus_v3\samples.tsv")
OUT_PATH = Path(r"data\discovery\taurus_v3\paired_for_run_discovery.tsv")

def main():
    # Read as whitespace-delimited with no header, then drop header row if present.
    df = pd.read_csv(IN_PATH, sep=r"\s+", header=None, engine="python")
    print("loaded", df.shape, "from", IN_PATH)

    if df.shape[1] < 6:
        raise ValueError(f"Expected at least 6 columns, found {df.shape[1]}")

    # If the first row looks like a header, drop it
    first_cell = str(df.iloc[0, 0]).strip().lower()
    if first_cell == "sample_id":
        df = df.iloc[1:].reset_index(drop=True)
        print("dropped header row detected in data")

    df = df.iloc[:, :6].copy()
    df.columns = ["sample_id", "subject_id", "disease", "site", "timepoint", "response"]

    out = df[["sample_id", "subject_id", "timepoint", "response"]].copy()
    out.to_csv(OUT_PATH, sep="\t", index=False)

    print("wrote", OUT_PATH)
    print(out.head(8).to_string(index=False))

if __name__ == "__main__":
    main()
