import warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
import dynamo as dyn
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from dynamo.preprocessing import Preprocessor
import scvelo as sc

adata = sc.read("object_velo.h5ad",cache=True)

celltype_key = "subtype"

preprocessor = Preprocessor()
preprocessor.preprocess_adata(adata, recipe="seurat")

dyn.tl.reduceDimension(adata)
dyn.pl.umap(adata, color='subtype')

dyn.tl.dynamics(adata, model='stochastic', est_method='negbin', cores=24)
dyn.tl.gene_wise_confidence(adata, group='subtype', lineage_dict={'MDSC_Chil3': ['macro_Mki67']})
dyn.tl.cell_velocities(adata, method='cosine', other_kernels_dict={'transform': 'sqrt'})
dyn.tl.cell_wise_confidence(adata)
dyn.tl.confident_cell_velocities(adata, group='subtype', lineage_dict={'MDSC_Chil3': ['macro_Mki67']})
dyn.pl.cell_wise_vectors(adata, color=['subtype'], inverse=False, quiver_length=1, quiver_size=1, pointsize=0.1, show_arrowed_spines=True)
dyn.pl.streamline_plot(adata, color=['subtype'],figsize=(5, 4),show_legend='on data',
                              color_key=["#824880","#cca6bf","#cc7eb1","#5b86ab",
                                         "#d5c666","#47885e","#d1edcb"])

dyn.vf.VectorField(adata, basis='umap', M=1000, pot_curl_div=True)
dyn.pl.plot_energy(adata, basis='umap')
dyn.pl.topography(adata, basis='umap', n=200, background='white', color=['subtype'], streamline_color='black', show_legend='on data', frontier=True)
dyn.pl.umap(adata,  color='umap_ddhodge_potential', frontier=True)

dyn.tl.cell_velocities(adata, basis='pca')
dyn.vf.VectorField(adata, basis='pca')
dyn.vf.speed(adata, basis='pca')
dyn.vf.curl(adata, basis='umap')
dyn.vf.divergence(adata, basis='pca')
dyn.vf.acceleration(adata, basis='pca')
dyn.vf.curvature(adata, basis='pca')

progenitor = adata.obs_names[adata.obs.subtype.isin(['MDSC_Chil3'])]
dyn.pd.fate(adata, basis='umap', init_cells=progenitor, interpolation_num=100,  direction='forward',
   inverse_transform=False, average=False, cores=24)

%%capture
fig, ax = plt.subplots()
ax = dyn.pl.topography(adata, color='subtype', ax=ax, save_show_or_return='return')

%%capture
instance = dyn.mv.StreamFuncAnim(adata=adata, color='subtype', ax=ax)

dyn.configuration.set_figure_params("dynamo", background="white")
dyn.pd.state_graph(adata, group='subtype', basis='pca', method='vf')
dyn.pl.streamline_plot(adata, color=['subtype'],figsize=(5, 4),show_legend='on data',
                              color_key=["#824880","#cca6bf","#cc7eb1","#5b86ab",
                                         "#d5c666","#47885e","#d1edcb"], pointsize=0.4,
                       save_show_or_return='save',save_kwargs={"path":"./", "prefix": 'Stream_tumor', "ext": 'svg', "dpi": None, "transparent": True, "close": True, "verbose": True})

dyn.pl.state_graph(adata,color=['subtype'],group='subtype',
                   color_key=["#824880","#cca6bf","#cc7eb1","#5b86ab",
                                         "#d5c666","#47885e","#d1edcb"],
                   basis='umap',figsize=[10,7],show_legend='on data',method='vf',
                  save_show_or_return='save',save_kwargs={"path":"./", "prefix": 'State_tumor', "ext": 'pdf', "dpi": None, "transparent": True, "close": True, "verbose": True})

fig, ax = plt.subplots(figsize=(12, 13))
sc.settings.set_figure_params("scvelo", vector_friendly = False)
dyn.pl.state_graph(adata,color=['subtype'],ax=ax,group='subtype',
                   color_key=["#824880","#cca6bf","#cc7eb1","#5b86ab",
                                         "#d5c666","#47885e","#d1edcb"],basis='umap',
                   show_legend=None,method='vf',edge_scale=0.5)
# layout
ax.set_title('')
plt.tight_layout()
fig.savefig(f'./State_tumor.png', bbox_inches='tight')

dyn.pl.state_graph(adata,color=['subtype'],group='subtype',
                   color_key=["#824880","#cca6bf","#cc7eb1","#5b86ab",
                                         "#d5c666","#47885e","#d1edcb"],basis='umap',
                   show_legend=None,method='vf')
