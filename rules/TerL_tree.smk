# Snakemake script
# 2.2.3 catch the terminase in the blastp fmt: self and neibours
rule find_terL:
    input: 
        blastp_fmt = "output/{sample}/blastp_fmt.tsv",
        faa = "output/{sample}/{sample}.faa"
    output: 
        terl_list = temp("output/{sample}/TerL.list"),
        terl_self = temp("output/{sample}/TerL_self.faa"),
        terl_neib = temp("output/{sample}/TerL_neib.faa"),
        terl_neib_fmt = temp("output/{sample}/TerL_neib_fmt.faa"),
        terl_faa = "output/{sample}/TerL.faa"
    params:
        totalname_dict = f'{db_path}/{db_prefix}_vConTACT2_proteins_totalname.pydict',
        nameseq_dict = f'{db_path}/{db_prefix}_vConTACT2_proteins_nameseq.pydict',
        genome_totalname_dict = f'{db_path}/{db_prefix}_genomes_totalname.pydict'
    log: f"{log_dir}/" + "{sample}_terL_tree.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''python {script_dir}/get_terL_self.py \\
    -i {input.blastp_fmt} -o {output.terl_list} \\
    -f {input.faa} -of {output.terl_self} -s {wildcards.sample} > {log}
python {script_dir}/get_seqs_from_dict.py \\
    -i {output.terl_list} -o {output.terl_neib} \\
    -ns {params.nameseq_dict} -tn {params.totalname_dict} >> {log}
python {script_dir}/terL_neib_id2name.py \\
    -i {output.terl_neib} -tn {params.genome_totalname_dict} \\
    -o {output.terl_neib_fmt} >> {log}
python {script_dir}/cat_terL.py \\
    --self {output.terl_self} --neib {output.terl_neib_fmt} \\
    -o {output.terl_faa} -s {wildcards.sample} >> {log}
'''


# 2.2.4 phylogenetic_tree
rule MAFFT_iqtree:
    input: "output/{sample}/TerL.faa"
    output:
        aln = "output/{sample}/TerL_tree/TerL.aln",
        tree = "output/{sample}/TerL_tree/TerL.aln.treefile"
    threads: 20
    params: wkdir = "output/{sample}/TerL_tree"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    log: f"{log_dir}/" + "{sample}_terL_tree.log"
    shell: '''if [ ! -d {params.wkdir} ];then mkdir -p {params.wkdir};fi
if [ $(grep -c '>' {input}) -gt 3 ];then 
    mafft --auto --quiet {input} > {output.aln}
    iqtree -s {output.aln} -T {threads} -B 1000 >> {log}
else
    echo "Error: {input} is not ok!"
    for i in {output};do
        touch ${{i}}
    done
fi
'''


# 2.2.5 ggtree visualization
rule phylotree_visualize:
    input: "output/{sample}/TerL_tree/TerL.aln.treefile"
    output: protected("output/{sample}/TerL.pdf")
    log: f"{log_dir}/" + "{sample}_terL_tree.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''
if [ -s {input} ];then
    python {script_dir}/plot_tree.py {input} {output} >> {log}
else
    touch {output}
fi
'''
