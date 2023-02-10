# Snakemake script

configfile: "config.yaml"
workdir: config['workdir']
db_path = config['db_path']
db_prefix = config['db_prefix']
script_dir = config['script_dir'].rstrip('\\/')

rule all:
    input:
        'Download_ok',
        "pydict_ok",
        "vc_accel_db"

rule download_db:
    output: temp('Download_ok')
    shell: '''cd {db_path}
function download(file){{
    Target={output.pfammap}
    success=0
    Ret=$(wget --tries 10 --retry-connrefused --waitretry=60 --timeout=60 -q http://inphared.s3.climb.ac.uk/${db_prefix}${{file}} || echo "404 Not Found")
    if echo $Ret | grep -i -q "404 Not Found"; then
        echo "server not responding for $file"
        exit 1
    else
        success=1
    fi
}}
for file in (
vConTACT2_proteins.faa
vConTACT2_family_annotations.tsv
vConTACT2_gene_to_genome.csv
vConTACT2_genus_annotations.tsv
vConTACT2_host_annotations.tsv
vConTACT2_lowest_taxa_annotations.tsv
vConTACT2_subfamily_annotations.tsv
genomes.fa
data.tsv
);do
    download ${{file}}
done
touch {output}
'''

rule make_py_dict:
    output: 
        temp(touch("pydict_ok"))
    shell: '''cd {db_path}
python scripts/mkdict.py {db_prefix}_genomes.fa
python scripts/mkdict.py {db_prefix}_vConTACT2_protein.fa
'''

rule vConTACT2_accelerate_db:
    output: temp(touch("vc_accel_db"))
    threads: 60
    shell: '''cd {db_path}
diamond makedb --in {db_prefix}_vConTACT2_proteins.faa \\
    -p {threads} --db {db_prefix}_vConTACT2_proteins.dmnd
diamond blastp --query {db_prefix}_vConTACT2_proteins.faa \\
    --db {db_prefix}_vConTACT2_proteins.dmnd -p {threads} \\
    --sensitive -o allVSall.dm.tsv
'''