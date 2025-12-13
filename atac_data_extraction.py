
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import pandas as pd

adata = sc.read_h5ad("GSM7876980_CEMBA190711_8E_rm_dlt.h5ad")
X = adata.X.T
X = X.astype(np.int32)
mmwrite("atac_matrix.mtx", X)

# Save features
pd.DataFrame(adata.var_names).to_csv("atac_features.tsv", sep="\t", index=False, header=False)

# Save barcodes
pd.DataFrame(adata.obs_names).to_csv("atac_barcodes.tsv", sep="\t", index=False, header=False)

# Save metadata
adata.obs.to_csv("atac_metadata.csv")
