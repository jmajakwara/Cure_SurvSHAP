#!/bin/sh
#SBATCH --partition=bigbatch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-23:00:00
#SBATCH --job-name=MC1000_Sim_log_low.R
#SBATCH --output=MC1000_Sim_log_low.txt
#SBATCH --chdir=/home-mscluster/jmajakwara/2024/SHAP/Simulation/Sample_Properties/	

R CMD BATCH --slave MC1000_Sim_log_low.R