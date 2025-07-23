# cpa-Predicting-combinatorial-drug-perturbations


## Conda Environment

open anaconda prompt on computer

create a python environment using the following command

```python
conda create -n cpa python=3.9
```

activate the environment using:
```python
conda activate cpa
```
deactivate environment
```python
conda deactivate cpa
```

## CPU only
```python
pip install torch==1.13.1+cpu --extra-index-url https://download.pytorch.org/whl/cpu
```

install latest version of CPA through our Github:

```python
pip install git+https://github.com/theislab/cpa

pip install scanpy
```


## for GPU
deactivate environment
```python
pip install torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117
```


## if the kernel is not shown in the Notebook run the following in the environment. 
Change the name-""
pip install ipykernel
python -m ipykernel install --user --name=cpa

## if the tkinter is not Found in 
import cpa
apt-get install python3.9-tk


Run this command in the repository folder to check if pretrained models exist:

ls -R | grep ".h5\|.h5ad\|.pt"

## After installation
Continue on a Jupiter notebook be sure to use the right kernel, choose it from the right. For subsampling use the scbasecount-py environment and for running optuna and cpa training (either case) use cpa kernel


## Subsampling of Tahoe-100M dataset

Create an environmnet to use for subsampling the Tahoe dataset.
You can use the conda environment provided in this git repository (CPA\environments\subsample-environment.yml).  
-Navigate to the folder where subsample-environment.yml is located.  
-Create the environment.
```python
conda env create -f subsample-environment.yml
```
Once the environment is created, activate it:
```
conda activate subsample-environment.yml
```
Activate the kernel of this environment to use for:
 [Tahoe_subsample_for_cell_pred.ipynb](https://github.com/mardiam/CPA/blob/main/Tahoe/Tahoe_subsample_for_cell_pred.ipynb) and [Tahoe_subsample_for_RDkit.ipynb](https://github.com/mardiam/CPA/blob/main/Tahoe/Tahoe_subsample_for_RDkit.ipynb)



 ## Optuna installation  
Run on the anaconda prompt, in the cpa environment (conda activate cpa) the following:  
 ```python 
 pip install optuna==2.10.1
 ```
# CPA Training

### Predicting perturbation responses for unseen cell-types (context transfer)


-Open [cpa_training_tahoe_on_cell_pred.ipynb](https://github.com/mardiam/CPA/blob/main/Tahoe/cpa_training_tahoe_on_cell_pred.ipynb)  
-Activate CPA kernel
-If you didn't subsample download on the working directory the [dataset](https://drive.google.com/file/d/1AH-ijqwa7rFyuxGeuv3ywyChhxXEKWby/view?usp=sharing).  
-Train the model or use the pretrained model you can download [from google drive]()


### Predicting drug perturbations using RDKit embeddings for drugs
-Open [cpa_training_tahoe_with_RDkit.ipynb](https://github.com/mardiam/CPA/blob/main/Tahoe/cpa_training_tahoe_%20with_RDkit.ipynb).

-Activate CPA kernel
-If you didn't subsample download on the working directory the [dataset](https://drive.google.com/file/d/18l3hwc9QKg0YrxRpESREFH0kCfNFFdKH/view?usp=sharing).   
-Train the model or use the pretrained [model](https://drive.google.com/drive/folders/1zjlzA-Js-FmRaSAX-6Jw5-q--F8ZxeXq?usp=sharing) you can download [from google drive](https://drive.google.com/drive/folders/1j12kn3rub9R981nxr_lRFwZC-kNz5t00?usp=sharing).
