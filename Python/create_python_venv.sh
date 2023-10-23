#!/bin/bash 

WD=/blue/rbacher/j.leary/py_envs  
cd $WD
module load conda
mamba create -p ${WD}/scLANE_sim_venv python=3.8
mamba activate scLANE_sim_venv
mamba install numpy 'pandas<2' 'matplotlib<3.7' scanpy scvelo
mamba deactivate
