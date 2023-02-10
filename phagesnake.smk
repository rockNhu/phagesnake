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

Samples, = glob_wildcards(os.path.join(fna_dir,"{name}.fasta"))
Samples = [sam for sam in Samples if '/' not in sam]


rule all:
    input:
        'Done-ANI',
        'Done-terL-tree',
        'Done-anno',
        'Done-vConTACT',
        'Done-Stat'


include: 'rules/nucl_align.smk'
include: 'rules/terL_tree.smk'
include: 'rules/annotations.smk'
include: 'rules/run_vConTACT2.smk'
include: 'rules/genome_stat.smk'

rule nucl_align:
    input: expand("output/{sample}/ANI_output",sample=Samples)
    output: temp('Done-ANI')
    shell: '''echo "Nucleotide alignment - finished. MMseqs2 + pyANI used."
touch {output}
'''

rule annotations:
    input:
        expand("output/{sample}/{sample}.png",sample=Samples),
        expand("output/{sample}/{sample}.gbk",sample=Samples)
    output: temp('Done-anno')
    shell: '''echo "TerL phylotree - finished."
echo "Prodigal + EggNog + DIAMOND + Biopython + dna_feature_viewer used."
touch {output}
'''

rule terL_tree:
    input: expand("output/{sample}/TerL.pdf",sample=Samples)
    output: temp('Done-terL-tree')
    shell: '''echo "TerL phylotree - finished. MAFFT + IQ-TREE + Biopython used."
touch {output}
'''

rule run_vConTACT:
    input:
        # expand("output/{sample}/c1.ntw",sample=Samples),
        # expand("output/{sample}/genome_by_genome_overview.csv",sample=Samples),
        expand("output/{sample}/{sample}_vConTACT.pdf",sample=Samples)
    output: temp('Done-vConTACT')
    shell: '''echo "vConTACT cluster - finished. vConTACT2 + graphanalyzer used."
touch {output}
'''

rule genome_stat:
    input: "seq_info.tsv"
    output: temp('Done-Stat')
    shell: '''echo "Genome statistic - finished. Only python3 used."
touch {output}
'''
onerror:
    sys.stderr.write('\n\n')
    sys.stderr.write(
        '*** An error occurred. Refer to the log file of the failed '
        'command for troubleshooting\n'
    )
    sys.stderr.write(
        'Issues can be raised at: https://github.com/rocknhu/phagesnake/issues\n'
    )