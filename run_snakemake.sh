source ~/.bashrc
conda activate phagesnake
nohup snakemake -s phagesnake.smk -c60 > s.log &
