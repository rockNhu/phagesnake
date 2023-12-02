# Snakemake script
# The nucleotide alignment sub-workflow

rule MMseqs_blastn:
    input:
        fa = fmt_fna_dir + "/{sample}.fasta",
        db = f'{db_path}/{db_prefix}_genomes.fa'
    output:
        tmp_dir = temp(directory('tmp/{sample}')),
        blastn_out = temp('output/{sample}/1.nucleotide_alignment/blastn.tsv')
    threads: 60
    log: f"{log_dir}/" + "{sample}_nucl_align.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 1.1 blastn to find the similar species genomes
echo "{wildcards.sample} : nucl_align start" > {log}
mmseqs easy-search --search-type 3 {input.fa} {input.db} {output.blastn_out} \\
    {output.tmp_dir} --threads {threads} --format-output "query,target,pident,qcov" >> {log}
'''

rule catch_nucl_neibours:
    input: 
        tsv = 'output/{sample}/1.nucleotide_alignment/blastn.tsv'
    output: 
        l = temp('output/{sample}/1.nucleotide_alignment/blastn.list'),
        fmt_blastn = 'output/{sample}/1.nucleotide_alignment/fmt_blastn.tsv'
    params:
        name_db = f'{db_path}/{db_prefix}_genomes_totalname.pydict',
        idt = 75,
        cov = 75
    log: f"{log_dir}/" + "{sample}_nucl_align.log"
    shell: '''# 1.2 parse the blastn output
python {script_dir}/get_neibour_nucl.py -i {input.tsv} \\
    -o_list {output.l} -o_fmt {output.fmt_blastn} -nucl_id_db {params.name_db} \\
    --identities {params.idt} --coverages {params.cov} >> {log}
'''

rule catch_neibours_fna:
    input: 
        fa = fmt_fna_dir + "/{sample}.fasta",
        ex_list = 'output/{sample}/1.nucleotide_alignment/blastn.list',
        totalname_dict = f'{db_path}/{db_prefix}_genomes_totalname.pydict',
        nameseq_dict = f'{db_path}/{db_prefix}_genomes_nameseq.pydict'
    output: 
        directory("output/{sample}/1.nucleotide_alignment/blastn_output")
    log: f"{log_dir}/" + "{sample}_nucl_align.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 1.3 catch the neibour fna and seperate neibours wtih format to "n_output"
if [ -s {input.ex_list} ];then
    python {script_dir}/get_seqs_from_dict.py -i {input.ex_list} \\
        -o {output} -ns {input.nameseq_dict} -tn {input.totalname_dict} --seperate >> {log}
    # add self fna
    cp {input.fa} {output}/{wildcards.sample}.fasta
else
    mkdir -p {output}
fi
'''

rule pyANI:
    input: 
        to_check = 'output/{sample}/1.nucleotide_alignment/blastn.list',
        nd = "output/{sample}/1.nucleotide_alignment/blastn_output"
    output:
        tab = "output/{sample}/1.nucleotide_alignment/pyani_out/ANIb_percentage_identity.tab"
    params:
        ani_dir = "output/{sample}/1.nucleotide_alignment/pyani_out",
        method = "ANIb"
    log: f"{log_dir}/" + "{sample}_nucl_align.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 1.4 do ANI using the pyani, the method is ANIb
if [ -s {input.to_check} ];then
    average_nucleotide_identity.py -i {input.nd} -o {params.ani_dir} -m {params.method} -g -f >> {log}
else
    echo "Error: empty neibours with {wildcards.sample}" >> {log}
    mkdir -p {params.ani_dir}
    touch {output.tab}
fi
if [ ! -f {output.tab} ];then touch {output.tab};fi
'''

rule ANI_plot:
    input:
        tab = "output/{sample}/1.nucleotide_alignment/pyani_out/ANIb_percentage_identity.tab"
    output:
        svg = "output/{sample}/1.nucleotide_alignment/ANIb_percentage_identity.svg",
        png = "output/{sample}/1.nucleotide_alignment/ANIb_percentage_identity.png"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell:'''# 1.5 replot the ANI heatmap
if [ -s {input.tab} ];then
    python {script_dir}/ani_heatmap.py -i {input.tab} -o {output.svg}
    python {script_dir}/ani_heatmap.py -i {input.tab} -o {output.png}
else
    touch {output.svg}
    touch {output.png}
fi
'''

rule nucl_align_all:
    input:
        svg = "output/{sample}/1.nucleotide_alignment/ANIb_percentage_identity.svg",
        png = "output/{sample}/1.nucleotide_alignment/ANIb_percentage_identity.png",
        blastn_out = "output/{sample}/1.nucleotide_alignment/blastn_output",
        fmt_blastn = "output/{sample}/1.nucleotide_alignment/fmt_blastn.tsv"
    output:
        temp(touch("output/{sample}_1.nucleotide_alignment"))