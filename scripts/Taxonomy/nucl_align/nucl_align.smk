# Snakemake script
# The nucleotide alignment sub-workflow
align_script_dir = script_dir + '/Taxonomy/nucl_align'

rule MMseqs_blastn:
    input:
        fa = fmt_fna_dir + "/{sample}.fasta",
        db = f'{db_path}/{db_prefix}_genomes.fa'
    output:
        tmp_dir = temp(directory('tmp/{sample}')),
        blastn_out = temp(homedir + '/output/{sample}/1.nucleotide_alignment/blastn.tsv')
    threads: 60
    log: f"{log_dir}/" + "{sample}_nucl_align.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 1 blastn to find the similar species genomes
echo "{wildcards.sample} : nucl_align start" > {log}
mmseqs easy-search --search-type 3 {input.fa} {input.db} {output.blastn_out} \\
    {output.tmp_dir} --threads {threads} --format-output "query,target,pident,qcov" >> {log}
'''

rule catch_nucl_neibours:
    input: 
        tsv = homedir + '/output/{sample}/1.nucleotide_alignment/blastn.tsv'
    output: 
        l = temp(homedir + '/output/{sample}/1.nucleotide_alignment/blastn.list'),
        fmt_blastn = homedir + '/output/{sample}/1.nucleotide_alignment/fmt_blastn.tsv'
    params:
        name_db = f'{db_path}/{db_prefix}_genomes_totalname.pydict',
        idt = 75,
        cov = 75
    log: f"{log_dir}/" + "{sample}_nucl_align.log"
    shell: '''# 2 parse the blastn output
python {align_script_dir}/get_neibour_nucl.py -i {input.tsv} \\
    -o_list {output.l} -o_fmt {output.fmt_blastn} -nucl_id_db {params.name_db} \\
    --identities {params.idt} --coverages {params.cov} >> {log}
'''

rule catch_neibours_fna:
    input: 
        fa = fmt_fna_dir + "/{sample}.fasta",
        ex_list = homedir + '/output/{sample}/1.nucleotide_alignment/blastn.list',
        totalname_dict = f'{db_path}/{db_prefix}_genomes_totalname.pydict',
        nameseq_dict = f'{db_path}/{db_prefix}_genomes_nameseq.pydict'
    output: 
        directory(homedir + "/output/{sample}/1.nucleotide_alignment/blastn_output")
    log: f"{log_dir}/" + "{sample}_nucl_align.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 3 catch the neibour fna and seperate neibours wtih format to "n_output"
if [ -s {input.ex_list} ];then
    python {align_script_dir}/get_seqs_from_dict.py -i {input.ex_list} \\
        -o {output} -ns {input.nameseq_dict} -tn {input.totalname_dict} --seperate >> {log}
    # add self fna
    cp {input.fa} {output}/{wildcards.sample}.fasta
else
    mkdir -p {output}
fi
'''

rule VIRIDIC_prepare:
    input: 
        to_check = homedir + '/output/{sample}/1.nucleotide_alignment/blastn.list',
        nd = homedir + "/output/{sample}/1.nucleotide_alignment/blastn_output"
    output:
        viridic_wkdir = directory(homedir + '/output/{sample}/1.nucleotide_alignment/viridic_out'),
        combined_neib = homedir + "/output/{sample}/1.nucleotide_alignment/viridic_out/combined_neibour.fasta"
    log: log_dir + "/{sample}_nucl_align.log"
    shell: '''# 4 catch all neibour genomes to one fasta file 
mkdir -p {output.viridic_wkdir}
if [ -s {input.to_check} ];then
    python {align_script_dir}/viridic_prepare.py -i {input.nd} -o {output.combined_neib}
else
    touch {output.combined_neib}
    echo "Error: empty neibours with {wildcards.sample}" >> {log}
fi
'''

rule VIRIDIC:
    input:
        viridic_wkdir = homedir + '/output/{sample}/1.nucleotide_alignment/viridic_out',
        fa = homedir + "/output/{sample}/1.nucleotide_alignment/viridic_out/combined_neibour.fasta"
    output:
        homedir + '/output/{sample}/1.nucleotide_alignment/viridic_out/done'
    params:
        viridic_script = homedir + '/external_software/viridic'
    log: log_dir + "/{sample}_nucl_align.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 5 alignment with other neibour genomes using VIRIDIC
if [ -s {input.fa} ];then
    cd {params.viridic_script}
    Rscript 00_viridic_master.R in={input.fa} projdir={input.viridic_wkdir} >> {log}
fi
touch {output}
'''

rule nucl_align_all:
    input:
        homedir + "/output/{sample}/1.nucleotide_alignment/viridic_out/done"
    output:
        temp(touch("output/{sample}_1.nucleotide_alignment"))
