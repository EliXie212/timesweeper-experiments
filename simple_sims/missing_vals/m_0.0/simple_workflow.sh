#!/bin/bash
#SBATCH --partition=dschridelab
#SBATCH --constraint=rhel8
#SBATCH --mem=8G
#SBATCH -c 16
#SBATCH --time=2:00:00
#SBATCH -J workflow
#SBATCH -o logfiles/workflow.%A.%a.out
#SBATCH -e logfiles/workflow.%A.%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lswhiteh@email.unc.edu
##SBATCH --array=0-2500:10
conda init bash
conda activate blinx
source activate blinx

srcdir=/proj/dschridelab/lswhiteh/timesweeper/timesweeper
configfile=config.yaml

#timesweeper condense --hft -o training_data.pkl yaml ${configfile}
#timesweeper train -i training_data.pkl --hft -n Missing_00 yaml ${configfile}
timesweeper plot_training -i training_data.pkl -n Missing_00 -o input_images