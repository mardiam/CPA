import cpa
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from collections import defaultdict
from tqdm import tqdm
import os

number_of_cells = 20000 # 5000, 10000, 15000, 20000

DATA_PATH = f'./resized_data_{number_of_cells}_cells.h5ad'
MODEL_PATH = f'./cpa_model{number_of_cells}/'
IMAGES_FOLDER = f'./images{number_of_cells}/'

# Create the images folder if it doesn't exist
if not os.path.exists(IMAGES_FOLDER):
    os.makedirs(IMAGES_FOLDER)

plt.ion()  # Enable interactive mode for matplotlib

print("Loading data...")
adata = sc.read(DATA_PATH)
adata.obs['split_1ct_MEC'].value_counts()
adata.X = adata.layers['counts'].copy()
print("Data loaded.")

print("Loading model...")
model = cpa.CPA.load(dir_path=MODEL_PATH, adata=adata, use_gpu=False)
print("Model loaded.")

#########################
### Trainning History ###
#########################
cpa.pl.plot_history(model)
# save to file
plt.savefig(IMAGES_FOLDER + 'history.png', dpi=300, bbox_inches='tight')
plt.show()
# plt.close('all')

#######################################
### Latent space UMAP visualization ###
#######################################
latent_outputs = model.get_latent_representation(adata, batch_size=1024)
sc.settings.verbosity = 3
latent_basal_adata = latent_outputs['latent_basal']
latent_adata = latent_outputs['latent_after']
sc.pp.neighbors(latent_basal_adata)
sc.tl.umap(latent_basal_adata)
latent_basal_adata
sc.pl.umap(latent_basal_adata, color=['condition_ID'], frameon=False, wspace=0.2)
plt.savefig(IMAGES_FOLDER + 'latent_basal_adata.png', dpi=300, bbox_inches='tight')
plt.show()
# plt.close('all')


# plt.close(fig='all')

sc.pp.neighbors(latent_adata)
sc.tl.umap(latent_adata)
sc.pl.umap(latent_adata, color=['condition_ID'], frameon=False, wspace=0.2)
plt.savefig(IMAGES_FOLDER + 'latent_adata.png', dpi=300, bbox_inches='tight')
plt.show()
# plt.close('all')

##################
### Evaluation ###
##################
model.predict(adata, batch_size=1024)
n_top_degs = [10, 20, 50, None] # None means all genes

results = defaultdict(list)
ctrl_adata = adata[adata.obs['condition_ID'] == 'CHEMBL504'].copy()
for cat in tqdm(adata.obs['cov_drug_dose'].unique()):
    if 'CHEMBL504' not in cat:
        cat_adata = adata[adata.obs['cov_drug_dose'] == cat].copy()

        deg_cat = f'{cat}'
        deg_list = adata.uns['rank_genes_groups_cov'][deg_cat]

        x_true = cat_adata.layers['counts'].toarray()
        x_pred = cat_adata.obsm['CPA_pred']
        x_ctrl = ctrl_adata.layers['counts'].toarray()

        x_true = np.log1p(x_true)
        x_pred = np.log1p(x_pred)
        x_ctrl = np.log1p(x_ctrl)

        for n_top_deg in n_top_degs:
            if n_top_deg is not None:
                degs = np.where(np.isin(adata.var_names, deg_list[:n_top_deg]))[0]
            else:
                degs = np.arange(adata.n_vars)
                n_top_deg = 'all'

            x_true_deg = x_true[:, degs]
            x_pred_deg = x_pred[:, degs]
            x_ctrl_deg = x_ctrl[:, degs]

            r2_mean_deg = r2_score(x_true_deg.mean(0), x_pred_deg.mean(0))
            r2_var_deg = r2_score(x_true_deg.var(0), x_pred_deg.var(0))

            r2_mean_lfc_deg = r2_score(x_true_deg.mean(0) - x_ctrl_deg.mean(0), x_pred_deg.mean(0) - x_ctrl_deg.mean(0))
            r2_var_lfc_deg = r2_score(x_true_deg.var(0) - x_ctrl_deg.var(0), x_pred_deg.var(0) - x_ctrl_deg.var(0))

            cov, cond, dose = cat.split('_')

            results['cell_type'].append(cov)
            results['condition'].append(cond)
            results['dose'].append(dose)
            results['n_top_deg'].append(n_top_deg)
            results['r2_mean_deg'].append(r2_mean_deg)
            results['r2_var_deg'].append(r2_var_deg)
            results['r2_mean_lfc_deg'].append(r2_mean_lfc_deg)
            results['r2_var_lfc_deg'].append(r2_var_lfc_deg)

df = pd.DataFrame(results)
df[df['n_top_deg'] == 20].to_csv(f'results{number_of_cells}.csv', index=False)
print(df[df['n_top_deg'] == 20])

####################################################
### We can further visualize these per condition ###
####################################################
for cat in adata.obs["cov_drug_dose"].unique():
    if "CHEMBL504" not in cat:
        cat_adata = adata[adata.obs["cov_drug_dose"] == cat].copy()

        cat_adata.X = np.log1p(cat_adata.layers["counts"].A)
        cat_adata.obsm["CPA_pred"] = np.log1p(cat_adata.obsm["CPA_pred"])

        deg_list = adata.uns["rank_genes_groups_cov"][f'{cat}'][:20]

        print(cat, f"{cat_adata.shape}")
        # cpa plot and save to file. do not show the plot
        cpa.pl.mean_plot(
            cat_adata,
            pred_obsm_key="CPA_pred",
            # path_to_save=None,
            deg_list=deg_list,
            # save_figure=True,
            show=False,
            verbose=True,
            path_to_save=IMAGES_FOLDER + f'{cat}_mean_plot.png'   
        )
        

######################################################
### Visualizing similarity between drug embeddings ###
######################################################
cpa_api = cpa.ComPertAPI(adata, model,
                         de_genes_uns_key='rank_genes_groups_cov',
                         pert_category_key='cov_drug_dose',
                         control_group='CHEMBL504',
                         )
cpa_plots = cpa.pl.CompertVisuals(cpa_api, fileprefix=None)
print(cpa_api.num_measured_points['train'])

drug_adata = cpa_api.get_pert_embeddings()
print(drug_adata.shape)

cpa_plots.plot_latent_embeddings(drug_adata.X, kind='perturbations', titlename='Drugs')
plt.savefig(IMAGES_FOLDER + 'drugs_latent_embeddings.png', dpi=300, bbox_inches='tight')
plt.show(block=True) 

