from pathlib import Path
import scanpy as sc
from anndata import AnnData


def load_ann_data(file_path: str | Path) -> AnnData:
    """
    Loads an AnnData object from a .h5ad or .csv file.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"No file found at {path}")

    if path.suffix == ".h5ad":
        return sc.read_h5ad(path)
    
    # Fallback to CSV loading
    return sc.read_csv(path)


def save_ann_data(adata: AnnData, file_path: str | Path) -> None:
    """
    Saves an AnnData object to a .h5ad file.
    """
    adata.write_h5ad(file_path)