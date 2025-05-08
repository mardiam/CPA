# cpa-Predicting-combinatorial-drug-perturbations


## Conda Environment

open anaconda prompt on computer

create a python environment using the following command

conda create -n cpa python=3.9

activate the environment using:

conda activate cpa

## CPU only

pip install torch==1.13.1+cpu --extra-index-url https://download.pytorch.org/whl/cpu

install latest version of CPA through our Github:

pip install git+https://github.com/theislab/cpa

pip install scanpy

## for GPU(κάρτα γραφικών)
pip install torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117

## for CPU(επεξεργαστής)
pip install torch==1.13.1+cpu --extra-index-url https://download.pytorch.org/whl/cpu
pip install git+https://github.com/theislab/cpa
pip install scanpy

## if the kernel is not shown in the Notebook run the following in the environment. 
Change the name-""
pip install ipykernel
python -m ipykernel install --user --name=cpa

## if the tkinter is not Found in 
import cpa
apt-get install python3.9-tk


Run this command in the repository folder to check if pretrained models exist:

ls -R | grep ".h5\|.h5ad\|.pt"


Continue on a Jupiter notebook that you open through  anaconda on desktop, be sure to use the right kernel, choose it from the right 

If you want to use the train.py (trains everysubsized dataset) and plot.py (you have to change the cell number one by one)
