source ~/.bashrc
conda activate phagesnake
nohup snakemake -s phagesnake.smk -k -c60 > s.log &
bash outzip.sh
