# Snakemake script
import os

configfile: "config.yaml"
db_path = config['db_path']
script_dir = os.path.abspath("scripts")
if not os.path.exists(db_path):
    os.makedirs(db_path)
workdir: db_path
db_prefix = config['db_prefix']

rule all:
    input:
        f"{db_prefix}_vConTACT2_proteins.faa",
        f"{db_prefix}_vConTACT2_family_annotations.tsv",
        f"{db_prefix}_vConTACT2_gene_to_genome.csv",
        f"{db_prefix}_vConTACT2_genus_annotations.tsv",
        f"{db_prefix}_vConTACT2_host_annotations.tsv",
        f"{db_prefix}_vConTACT2_lowest_taxa_annotations.tsv",
        f"{db_prefix}_vConTACT2_subfamily_annotations.tsv",
        f"{db_prefix}_genomes.fa",
        f"{db_prefix}_data.tsv",
        'allVSall.dm.tsv',
        f'{db_prefix}_genomes_totalname.pydict',
        f'{db_prefix}_genomes_nameseq.pydict',
        f'{db_prefix}_vConTACT2_proteins_totalname.pydict',
        f'{db_prefix}_vConTACT2_proteins_nameseq.pydict'


rule download_db:
    params:
        main_url = 'https://millardlab-inphared.s3.climb.ac.uk'
    output: 
        f"{db_prefix}_vConTACT2_proteins.faa",
        f"{db_prefix}_vConTACT2_family_annotations.tsv",
        f"{db_prefix}_vConTACT2_gene_to_genome.csv",
        f"{db_prefix}_vConTACT2_genus_annotations.tsv",
        f"{db_prefix}_vConTACT2_host_annotations.tsv",
        f"{db_prefix}_vConTACT2_lowest_taxa_annotations.tsv",
        f"{db_prefix}_vConTACT2_subfamily_annotations.tsv",
        f"{db_prefix}_genomes.fa",
        f"{db_prefix}_data.tsv"
    shell: '''
function download() {{
    file=$1
    if [ ! -s ${{file}}.gz ] && [ ! -s ${{file}} ];then
        echo "Downloading $file"
        Ret=$(wget --tries 5 --timeout=60 -q {params.main_url}/${{file}} || \\
            wget --tries 5 --timeout=60 -q {params.main_url}/${{file}}.gz || \\
            echo "404 Not Found")
        if echo $Ret | grep -i -q "404 Not Found"; then
            echo "server not responding for $file"
            exit 1
        else
            echo "$file was download successfully"
        fi
    else
        echo "$file was download successfully"
    fi
}}
for file in {output};do
    download ${{file}}
done
gunzip *.gz
'''

rule make_py_dict:
    input:
        genomes = f'{db_prefix}_genomes.fa',
        protein = f'{db_prefix}_vConTACT2_proteins.faa'
    output: 
        f'{db_prefix}_genomes_totalname.pydict',
        f'{db_prefix}_genomes_nameseq.pydict',
        f'{db_prefix}_vConTACT2_proteins_totalname.pydict',
        f'{db_prefix}_vConTACT2_proteins_nameseq.pydict',
    shell: '''
python {script_dir}/mkdict.py {input.genomes}
python {script_dir}/mkdict.py {input.protein}
'''

rule vConTACT2_accelerate_db:
    input:
        protein = f'{db_prefix}_vConTACT2_proteins.faa',
    output:
        dmnd = f'{db_prefix}_vConTACT2_proteins.dmnd',
        ava = 'allVSall.dm.tsv'
    threads: 60
    shell: '''
diamond makedb --in {input.protein} -p {threads} --db {output.dmnd}
diamond blastp --query {input.protein} --db {output.dmnd} -p {threads} \\
    --sensitive -o {output.ava}
'''

rule make_eggnog_database:
    output:
        f'{db_path}eggnog.db'
    shell: '''download_eggnog_data.py --data_dir {db_path} -M'''
