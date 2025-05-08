import cpa
import scanpy as sc
import pandas as pd
import numpy as np

def resize_dataset(data, n_cells):
    """Resizes the dataset to a target number of cells and saves it."""
    
    print("initial size: " + str(data.n_obs))

    n_cells = min(n_cells, data.n_obs)

    # Creates a random number generator (RNG) with a fixed seed (42) for reproducibility.
    # Randomly selects n_cells cell names from the dataset.
    # Create a copy of the new data, which were randomly selected, without modifying the original dataset.
    data = data[np.random.default_rng(seed=42).choice(data.obs_names, size=n_cells, replace=False)].copy()  

    # scanpy allows to store multiple versions of data, counts refers to raw data, which are used in CPA
    data.X = data.layers['counts'].copy()

    print("size after: " + str(data.n_obs))

    # Save the resized dataset
    save_path = f"./resized_data_{n_cells}_cells.h5ad"
    data.write(save_path)
    print(f"Resized dataset saved at: {save_path}")

    return data
    
    
def train(adata, save_path):
    """Train the CPA model."""
    
    adata.obs['split_1ct_MEC'].value_counts()
    adata.X = adata.layers['counts'].copy()
    
    cpa.CPA.setup_anndata(adata,
                        perturbation_key='condition_ID',
                        dosage_key='log_dose',
                        control_group='CHEMBL504',
                        batch_key=None,
                        is_count_data=True,
                        categorical_covariate_keys=['cell_type'],
                        deg_uns_key='rank_genes_groups_cov',
                        deg_uns_cat_key='cov_drug_dose',
                        max_comb_len=2,
                        )

    ae_hparams = {
        "n_latent": 128,
        "recon_loss": "nb",
        "doser_type": "logsigm",
        "n_hidden_encoder": 512,
        "n_layers_encoder": 3,
        "n_hidden_decoder": 512,
        "n_layers_decoder": 3,
        "use_batch_norm_encoder": True,
        "use_layer_norm_encoder": False,
        "use_batch_norm_decoder": True,
        "use_layer_norm_decoder": False,
        "dropout_rate_encoder": 0.1,
        "dropout_rate_decoder": 0.1,
        "variational": False,
        "seed": 434,
    }

    trainer_params = {
        "n_epochs_kl_warmup": None,
        "n_epochs_pretrain_ae": 30,
        "n_epochs_adv_warmup": 50,
        "n_epochs_mixup_warmup": 3,
        "mixup_alpha": 0.1,
        "adv_steps": 2,
        "n_hidden_adv": 64,
        "n_layers_adv": 2,
        "use_batch_norm_adv": True,
        "use_layer_norm_adv": False,
        "dropout_rate_adv": 0.3,
        "reg_adv": 20.0,
        "pen_adv": 20.0,
        "lr": 0.0003,
        "wd": 4e-07,
        "adv_lr": 0.0003,
        "adv_wd": 4e-07,
        "adv_loss": "cce",
        "doser_lr": 0.0003,
        "doser_wd": 4e-07,
        "do_clip_grad": False,
        "gradient_clip_value": 1.0,
        "step_size_lr": 45,
    }

    model = cpa.CPA(adata=adata,
                    split_key='split_1ct_MEC',
                    train_split='train',
                    valid_split='valid',
                    test_split='ood',
                    **ae_hparams,
                )

    model.train(max_epochs=500,
                use_gpu=False,
                batch_size=64,
                plan_kwargs=trainer_params,
                early_stopping_patience=30,
                check_val_every_n_epoch=5,
                save_path=save_path,
            )
    
    return model

def main():
    sc.settings.set_figure_params(dpi=100)

    data_path = './combo_sciplex_prep_hvg_filtered.h5ad'

    adata = sc.read(data_path)
    
    data_sizes = [5000, 10000, 15000, 20000]
    for size in data_sizes:
        print("----------------- Training model with " + str(size) + " cells --------------")
        resized_data = resize_dataset(adata, n_cells=size)
        model = train(resized_data, save_path='./cpa_model' + str(size))

        print("Model trained and saved for " + str(size) + " cells.")
        print("--------------------------------------------------")
main()
    

    