#Importing required libraries
import pandas as pd
from pathlib import Path
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
import scanpy as sc
from scipy import sparse
import scipy.io
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# Setting up cache as suggested by ats atlas (downloads metadata/expression on first run)
download_base = Path('./abc_data') 
abc_cache = AbcProjectCache.from_cache_dir(download_base)
print(abc_cache.current_manifest) 

# Loaing all the metadata with cluster annotation, there is an option to download without cluster annotaion with just cell_metadata
cell_extended = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='cell_metadata_with_cluster_annotation')
cell_extended.set_index('cell_label', inplace=True)

#Preparing and downloading the expression matrix from hippocampus through available directory and given label
#This data is big 7.42GB, we will be subsettiing this whole matrix for few specific cell types
rna_matrix_label = 'WMB-10Xv3-HPF'
rna_file = abc_cache.get_data_path(directory='WMB-10Xv3', file_name=f'{rna_matrix_label}/raw')
adata = sc.read_h5ad(rna_file, backed='r') 
print("Number of genes = ", len(adata.var))
adata.var.index[0:5]
print("Number of cells = ", len(adata.obs))
adata.obs.index[0:5]

#Identifying the available region in the downloaded cluster annotated cell metadata fromw where we will find acronym for different region of brain
cell_extended["region_of_interest_acronym"].unique()

# Filter for hippocampus data (HIP is used as acronym) and creating a separate copy to not modify raw downloaded data
hip_meta = cell_extended[cell_extended["region_of_interest_acronym"] == "HIP"].copy()

# Some exploratory data analysis
class_counts = hip_meta["class"].value_counts().reset_index()
class_counts.columns = ['Cell Class', 'Count']

# some depth exploratory data analysis for hippocampus cell types with dot plot
sns.set_style("white")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
fig, ax = plt.subplots(figsize=(6.5,7))
cmap = plt.cm.viridis
norm = Normalize(vmin=0, vmax=len(class_counts)-1)
colors = [cmap(norm(i)) for i in range(len(class_counts))]
sizes = (class_counts['Count'] / class_counts['Count'].max()) * 500 + 50
scatter = ax.scatter(
    x=class_counts['Count'],
    y=range(len(class_counts)),
    s=sizes,
    c=range(len(class_counts)),
    cmap='viridis',
    edgecolor='black',
    linewidth=1.5,
    alpha=0.8
)
ax.set_yticks(range(len(class_counts)))
ax.set_yticklabels(class_counts['Cell Class'], fontsize=8)
ax.tick_params(axis='x', labelsize=8)
ax.set_xlabel('Number of Cells', fontsize=12, fontweight='bold')
ax.set_ylabel('Cell Class', fontsize=12, fontweight='bold')
ax.set_title('Cell Class Distribution in Hippocampus', 
             fontsize=14, fontweight='bold', pad=15)
ax.grid(False)
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
x_margin = class_counts['Count'].max() * 0.05
ax.set_xlim(-x_margin, class_counts['Count'].max() + x_margin)
ax.invert_yaxis()
cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
cbar.set_label('Cell Class Index', fontsize=10, fontweight='bold')
cbar.ax.tick_params(labelsize=8)
#saving result as a pdf file
plt.savefig('hippocampus_cell_distribution.pdf', bbox_inches='tight')
plt.tight_layout()
plt.show()


# From the fig above we identifies different class/cell types in hippocampus, I decided to choose 04 DG-IMN Glut and 06 CTX-CGE GABA to subset for two cell types
processed_dir = Path("data/processed")
selected_classes = ["04 DG-IMN Glut", "06 CTX-CGE GABA"]
hip_meta_sub = hip_meta[hip_meta["class"].isin(selected_classes)].copy()
hip_meta_sub_reset = hip_meta_sub.reset_index() 

# saving output to csv for future use, so that we do not need to download the whole matrix and metadata again
meta_out = processed_dir / "hippocampus_selected_classes_metadata.csv"
processed_dir.mkdir(parents=True, exist_ok=True)
hip_meta_sub_reset.to_csv(meta_out, index=False)
print(hip_meta_sub["class"].value_counts())
print("Selected cells:", hip_meta_sub.shape[0])


#Matching metadata to adata.obs. cell_extended index is cell_label; adata.obs_names are also cell_label-like
selected_labels = hip_meta_sub.index.astype(str)
obs_index = adata.obs_names.astype(str)

mask = np.isin(obs_index, selected_labels)
print("Total HPF cells:", adata.n_obs)
print("Matching selected cells:", mask.sum())

# 3. Pull only those cells into memory from ABC suggestion
X_sub = adata.X[mask, :]
if hasattr(X_sub, "toarray"):
    X_sub = X_sub.toarray()
print("Subset matrix shape:", X_sub.shape)

# Extract gene_names and If from var
gene_names = adata.var["gene_symbol"].astype(str) 
print(gene_names.head()) 

# sparse genes x cells, Transposing to convert genes in row and cells in column
expr_sparse = sparse.csr_matrix(X_sub.T) 
processed_dir = Path("data/processed")
processed_dir.mkdir(parents=True, exist_ok=True)

#creating path to save the output required to create seurat object
mm_path = processed_dir / "hippocampus_expr.mtx"
genes_path = processed_dir / "genes.tsv"
barcodes_path = processed_dir / "barcodes.tsv"

#Saving result to specific directory
scipy.io.mmwrite(str(mm_path), expr_sparse)
pd.DataFrame(gene_names).to_csv(genes_path, sep="\t", header=False, index=False)

selected_labels_in_order = obs_index[mask]
pd.DataFrame(selected_labels_in_order).to_csv(
    barcodes_path, sep="\t", header=False, index=False
)
