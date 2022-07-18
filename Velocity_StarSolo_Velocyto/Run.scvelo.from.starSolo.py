#!/bin/python3
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import sparse
import anndata
import scanpy as sc
import scvelo as scv
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

path = sys.argv[1]
metadata = sys.argv[2]

X = sc.read_mtx(path+'Gene/raw/matrix.mtx')
X = X.X.transpose()
mtxU = np.loadtxt(path+'Velocyto/raw/unspliced.mtx', skiprows=3, delimiter=' ')
mtxS = np.loadtxt(path+'Velocyto/raw/spliced.mtx', skiprows=3, delimiter=' ')
mtxA = np.loadtxt(path+'Velocyto/raw/ambiguous.mtx', skiprows=3, delimiter=' ')
shapeU = np.loadtxt(path+'Velocyto/raw/unspliced.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
shapeS = np.loadtxt(path+'Velocyto/raw/spliced.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
shapeA = np.loadtxt(path+'Velocyto/raw/ambiguous.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
spliced = sparse.csr_matrix((mtxS[:,2], (mtxS[:,0]-1, mtxS[:,1]-1)), shape = shapeS).transpose()
unspliced = sparse.csr_matrix((mtxU[:,2], (mtxU[:,0]-1, mtxU[:,1]-1)), shape = shapeU).transpose()
ambiguous = sparse.csr_matrix((mtxA[:,2], (mtxA[:,0]-1, mtxA[:,1]-1)), shape = shapeA).transpose()
obs = pd.read_csv(path+'Velocyto/raw/barcodes.tsv', header = None, index_col = 0)
obs.index.name = None
var = pd.read_csv(path+'Velocyto/raw/features.tsv', sep='\t', names = ('gene_ids', 'feature_types'), index_col = 1)
adata = anndata.AnnData(X = X, obs = obs, var = var, layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous})
adata.var_names_make_unique()
selected_barcodes = pd.read_csv(path+'Gene/filtered/barcodes.tsv', header = None)
cell_meta = pd.read_csv(metadata)
adata = adata[cell_meta['cell']]
adata.obs = cell_meta
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.tl.umap(adata)
scv.tl.recover_dynamics(adata)
scv.tl.latent_time(adata)
adata.write(path+'scvelo.h5ad', compression = 'gzip')
scv.pl.proportions(adata,save=path+'proportion.png',show=None, groupby='PubName')
scv.pl.velocity_embedding_stream(adata, basis="X_umap", color="PubName",save=path+'UMAP_velocity.png',show=None)
scv.pl.velocity_embedding(adata, basis="X_umap", arrow_length=3, arrow_size=2, dpi=120, color="PubName",save=path+'UMAP_velocity.moreArrow.png',show=None)
scv.pl.scatter(adata, color="latent_time", color_map="gnuplot",save=path+'latent_time.png',show=None)
top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", n_convolve=100,save=path+'heatmap.png',show=None)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False,save=path+'Top.Gene.scatterplot.png',show=None)

