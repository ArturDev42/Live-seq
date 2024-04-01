import os
import pandas as pd
import scanpy as sc
import requests
from pathlib import Path

# Set root directory based on project file presence
root_dir = Path.cwd()  # Adjust this to your project's root directory logic

# Download count matrix if it doesn't exist
file_path = root_dir / "data/GSE141064_count.final.csv.gz"
if not file_path.exists():
    url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141064&format=file&file=GSE141064%5Fcount%2Efinal%2Ecsv%2Egz"
    r = requests.get(url)
    with open(file_path, 'wb') as f:
        f.write(r.content)

# Load data
count_all = pd.read_csv(file_path, index_col=0)
print(count_all.shape)
print(count_all)

meta_all = pd.read_csv(root_dir / "data/meta.final.csv", index_col=0)
meta_all.index = meta_all['sample_ID']
print(meta_all)
print(meta_all.columns)
print(meta_all.Extraction_time_h)
non_nan_rows_exraction_time = meta_all[meta_all['Extraction_time_h'].notnull()]
print(non_nan_rows_exraction_time)
print(non_nan_rows_exraction_time.Extraction_time_h)

# # Ensure matching sample_IDs
# # assert count_all.columns.equals(meta_all.index)

# # Create an AnnData object similar to Seurat object
# adata = sc.AnnData(count_all.T)
# adata.obs = meta_all

# # Add feature meta
# gene_name = pd.read_csv(
#     root_dir / "data/mouseGeneTable87_mCherry_EGFP.txt", sep="\t", index_col=0)
# gene_name_uni = gene_name[~gene_name.index.duplicated()]
# # assert all(adata.var_names.isin(gene_name_uni.index))
# adata.var['gene_ids'] = gene_name_uni.loc[adata.var_names, 'ensembl_gene_id']

# # Quality control plots can be created using scanpy's plotting functions, similar to Seurat's VlnPlot and FeatureScatter

# # Saving the processed object
# adata.write(root_dir / "01_preprocessing/Seu.all.h5ad")
