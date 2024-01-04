# snakemake script
#-* coding = UTF-8 *-
# @Author = Shixuan Huang (Rock Nhu)
import os
import datetime

def fmt_dir_end(d):
    return d.rstrip('\\/')


configfile: "config.yaml"
workdir: fmt_dir_end(config['workdir'])
log_dir = fmt_dir_end(config['log_dir'])
if not os.path.exists(log_dir):
    os.makedirs(log_dir)
Conda_env_dir = fmt_dir_end(config['Conda_yaml_dir'])
db_path = fmt_dir_end(config['db_path'])
db_prefix = config['db_prefix']
run_vc = config['run_vConTACT']
script_dir = os.path.join(fmt_dir_end(config['workdir']), 'scripts')
# fmt input sample files in fasta fmt
fna_dir = fmt_dir_end(config['fna_dir'])
fmt_fna_dir = fmt_dir_end(config['fmt_fna_dir'])
os.system(f'python {script_dir}/fmt_input_file.py -i {fna_dir} -o {fmt_fna_dir} > {log_dir}/fmt.log')
# the fmt input is gate of protocol
Samples, = glob_wildcards(os.path.join(fmt_fna_dir,"{name}.fasta"))
Samples = [sam for sam in Samples if '/' not in sam]
# current time to make statistic
start_time = datetime.datetime.now().strftime('%Y%m%d')

if run_vc == True:
    rule all:
        input:
            'Done-all'

    include: 'rules/run_vConTACT2.smk'
else:
    rule all:
        input:
            'Done-all-none-vcontact'

include: 'rules/nucl_align.smk'
include: 'rules/TerL_tree.smk'
include: 'rules/annotations.smk'
include: 'rules/genome_stat.smk'


rule nucl_align:
    input:
        expand("output/{sample}_1.nucleotide_alignment",sample=Samples),
    output: temp('Done-ANI')
    shell: '''echo "Nucleotide alignment - finished. MMseqs2 + pyANI used."
touch {output}
'''

rule annotations:
    input: 
        expand("output/{sample}_2.annotations-clean",sample=Samples),
        expand("output/{sample}_3.secuity_check-clean",sample=Samples)
    output: temp('Done-anno')
    shell: '''echo "Annotations - finished."
echo "Prodigal + EggNog + DIAMOND + Biopython + dna_feature_viewer used."
touch {output}
'''

rule TerL_tree:
    input: 
        expand("output/{sample}_4.TerL_phylogenetic_tree",sample=Samples),
    output: temp('Done-TerL-tree')
    shell: '''echo "TerL phylotree - finished. MAFFT + IQ-TREE + Biopython used."
touch {output}
'''

rule genome_stat: 
    input: f"output/5.seq_info{start_time}.tsv"
    output: temp('Done-Stat')
    shell: '''echo "Genome statistic - finished. Only python3 used."
touch {output}
'''

rule run_vConTACT:
    input: 
        expand('output/{sample}_6.vConTACT2_network', sample=Samples)
    output: temp('Done-vConTACT')
    shell: '''echo "vConTACT cluster - finished. vConTACT2 + graphanalyzer used."
touch {output}
'''

rule finished_clean:
    input:
        'Done-ANI',
        'Done-TerL-tree',
        'Done-anno',
        'Done-Stat',
        'Done-vConTACT'
    output:
        temp('Done-all')
    shell: '''# move the finished input out
if [ ! -d finished_fna ];then mkdir finished_fna;fi
if [ ! -d finished_fmt_fna ];then mkdir finished_fmt_fna;fi
if [ $(ls -A {fna_dir}) ];then
    mv {fna_dir}/* finished_fna
    mv {fmt_fna_dir}/* finished_fmt_fna
fi
touch {output}
'''

rule no_vc_finished_clean:
    input:
        'Done-ANI',
        'Done-TerL-tree',
        'Done-anno',
        'Done-Stat'
    output:
        temp('Done-all-none-vcontact')
    shell: '''# move the finished input out
if [ ! -d finished_fna ];then mkdir finished_fna;fi
if [ ! -d finished_fmt_fna ];then mkdir finished_fmt_fna;fi
if [ $(ls -A {fna_dir}) ];then
    mv {fna_dir}/* finished_fna
    mv {fmt_fna_dir}/* finished_fmt_fna
fi
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
