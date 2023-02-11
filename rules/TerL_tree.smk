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
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''python {script_dir}/get_terL_self.py \\
    -i {input.blastp_fmt} -o {output.terl_list} \\
    -f {input.faa} -of {output.terl_self} -s {wildcards.sample}
python {script_dir}/get_seqs_from_dict.py \\
    -i {output.terl_list} -o {output.terl_neib} \\
    -ns {params.nameseq_dict} -tn {params.totalname_dict}
python {script_dir}/terL_neib_id2name.py \\
    -i {output.terl_neib} -tn {params.genome_totalname_dict} \\
    -o {output.terl_neib_fmt}
python {script_dir}/cat_terL.py \\
    --self {output.terl_self} --neib {output.terl_neib_fmt} \\
    -o {output.terl_faa} -s {wildcards.sample}
'''


# 2.2.4 phylogenetic_tree
rule MAFFT_iqtree:
    input: "output/{sample}/TerL.faa"
    output: "output/{sample}/TerL.aln", "output/{sample}/TerL.aln.treefile"
    threads: 20
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''if [ $(grep -c '>' {input}) -gt 3 ];then 
    mafft --auto {input} > {output[0]}
    iqtree -s {output[0]} -T {threads} -B 1000
else
    echo "Error: {input} is not ok!"
    touch {output[0]}
    touch {output[1]}
fi
'''


# 2.2.5 ggtree visualization
rule phylotree_visualize:
    input: "output/{sample}/TerL.aln.treefile"
    output: "output/{sample}/TerL.pdf"
    params:
        wkdir = "output/{sample}"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''
if [ -s {input} ];then
    python {script_dir}/plot_tree.py {input} {output}
    # cleanup
    if [ ! -d {params.wkdir}/temp ];then mkdir {params.wkdir}/temp;fi
    mv {params.wkdir}/TerL.aln* {params.wkdir}/temp
else
    touch {output}
fi
'''
