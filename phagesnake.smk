# snakemake script
#-* coding = UTF-8 *-
# @Author = Shixuan Huang (Rock Nhu)
import os


configfile: "config.yaml"
workdir: config['workdir']
fna_dir = config['fna_dir'].rstrip('\\/')
db_path = config['db_path']
db_prefix = config['db_prefix']
script_dir = config['script_dir'].rstrip('\\/')
py_dict_dir = config['py_dict_dir']

Samples, = glob_wildcards(os.path.join(fna_dir,"{name}.fasta"))
Samples = [sam for sam in Samples if '/' not in sam]

rule all:
    input:
        expand("output/{sample}/ANI_output",sample=Samples),
        expand("output/{sample}/{sample}.png",sample=Samples),
        expand("output/{sample}/TerL.pdf",sample=Samples),
        expand("output/{sample}/{sample}.gbk",sample=Samples),
        expand("output/{sample}/c1.ntw",sample=Samples),
        expand("output/{sample}/genome_by_genome_overview.csv",sample=Samples),
        #expand("output/{sample}/{sample}_vConTACT.pdf",sample=Samples),
        "seq_info.tsv"


include: 'rules/nucl_align.smk'
include: 'rules/terL_tree.smk'
include: 'rules/annotations.smk'
include: 'rules/genome_stat.smk'