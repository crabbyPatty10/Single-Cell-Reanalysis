import anndata as ad

adata = ad.read_csv(
    "data.csv",          # i dont think we have the dataset in right now, but this is how you would load it
    delimiter=",",
    first_column_names=None,  # set to True if first col is obs names
    dtype="float32"
)

def load_data(file_path):
    with open(file_path, 'r') as file:
        data = file.read()
    return data