# Snakemake script
rule download_db:
    output: temp('Download_ok')
    shell: '''cd {db_path}
wget http://inphared.s3.climb.ac.uk/${db_prefix}_vConTACT2_proteins.faa
wget http://inphared.s3.climb.ac.uk/${db_prefix}_vConTACT2_family_annotations.tsv
wget http://inphared.s3.climb.ac.uk/${db_prefix}_vConTACT2_gene_to_genome.csv
wget http://inphared.s3.climb.ac.uk/${db_prefix}_vConTACT2_genus_annotations.tsv
wget http://inphared.s3.climb.ac.uk/${db_prefix}_vConTACT2_host_annotations.tsv
wget http://inphared.s3.climb.ac.uk/${db_prefix}_vConTACT2_lowest_taxa_annotations.tsv
wget http://inphared.s3.climb.ac.uk/${db_prefix}_vConTACT2_subfamily_annotations.tsv
wget http://inphared.s3.climb.ac.uk/${db_prefix}_genomes.fa
wget http://inphared.s3.climb.ac.uk/${db_prefix}_data.tsv
touch {output}
'''

rule make_py_dict:
    output: 
        "1Oct2022_vConTACT2_proteins_totalname.pydict"
    shell: '''cd {db_path}
python scripts/mkdict.py {db_prefix}_genomes.fa
python scripts/mkdict.py {db_prefix}_vConTACT2_protein.fa
'''

rule vConTACT2_accelerate_db:
    output: "allVSall.dm.tsv"
    threads: 60
    shell: '''cd {db_path}
diamond makedb --in {db_prefix}_vConTACT2_proteins.faa \\
    -p {threads} --db {db_prefix}_vConTACT2_proteins.dmnd
diamond blastp --query {db_prefix}_vConTACT2_proteins.faa \\
    --db {db_prefix}_vConTACT2_proteins.dmnd -p {threads} \\
    --sensitive -o {output}
'''