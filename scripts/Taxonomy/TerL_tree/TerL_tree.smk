# Snakemake script
# The phylogenetic tree based on TerL sub-workflow
TerL_script_dir = script_dir + "/Taxonomy/TerL_tree"
tree_based_on = config[tree_based_on]

rule find_terL:
    input: 
        blastp_fmt = homedir + "/output/{sample}/1.annotations/blastp_fmt.tsv",
        faa = homedir + "/output/{sample}/1.annotations/{sample}.faa"
    output: 
        terl_list = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL_tree/TerL.list",
        terl_self = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL_tree/TerL_self.faa",
        terl_neib = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL_tree/TerL_neib.faa",
        terl_neib_fmt = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL_tree/TerL_neib_fmt.faa",
        terl_faa = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL_tree/TerL.faa"
    params:
        totalname_dict = f'{db_path}/{db_prefix}_vConTACT2_proteins_totalname.pydict',
        nameseq_dict = f'{db_path}/{db_prefix}_vConTACT2_proteins_nameseq.pydict',
        genome_totalname_dict = f'{db_path}/{db_prefix}_genomes_totalname.pydict'
    log: log_dir + "/{sample}_terL_tree.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 1 catch the terminase in the blastp fmt: self and neibours
echo "{wildcards.sample} : TerL_tree start" > {log}
python {TerL_script_dir}/get_terL_self.py \\
    -i {input.blastp_fmt} -o {output.terl_list} \\
    -f {input.faa} -of {output.terl_self} -s {wildcards.sample} {tree_based_on} >> {log}
python {script_dir}/get_seqs_from_dict.py \\
    -i {output.terl_list} -o {output.terl_neib} \\
    -ns {params.nameseq_dict} -tn {params.totalname_dict} >> {log}
python {TerL_script_dir}/terL_neib_id2name.py \\
    -i {output.terl_neib} -tn {params.genome_totalname_dict} \\
    -o {output.terl_neib_fmt} >> {log}
python {TerL_script_dir}/cat_terL.py \\
    --self {output.terl_self} --neib {output.terl_neib_fmt} \\
    -o {output.terl_faa} -s {wildcards.sample} >> {log}
'''

rule MAFFT_iqtree:
    input: homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL_tree/TerL.faa"
    output:
        aln = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL_tree/TerL.aln",
        tree = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL.aln.treefile"
    threads: 20
    params: wkdir = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL_tree"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    log: log_dir + "/{sample}_terL_tree.log"
    shell: '''# 2 phylogenetic_tree
if [ ! -d {params.wkdir} ];then mkdir -p {params.wkdir};fi
if [ $(grep -c '>' {input}) -gt 3 ];then 
    mafft --auto --quiet {input} > {output.aln}
    iqtree -s {output.aln} -T {threads} -B 1000 --redo-tree >> {log}
    cp {params.wkdir}/TerL.aln.treefile {output.tree}
else
    echo "Error: {input} is not ok!"
    for i in {output};do
        touch ${{i}}
    done
fi
'''

rule phylotree_visualize:
    input: homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL.aln.treefile"
    output:
        svg_out = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL.svg",
        png_out = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL.png"
    log: log_dir + "/{sample}_terL_tree.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 3 ggtree visualization
if [ -s {input} ];then
    python {TerL_script_dir}/plot_tree.py {input} {output.svg_out} >> {log}
    python {TerL_script_dir}/plot_tree.py {input} {output.png_out} >> {log}
else
    touch {output}
fi
'''

rule TerL_clean:
    input:
        terl_tree = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL.aln.treefile",
        svg_out = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL.svg",
        png_out = homedir + "/output/{sample}/4.TerL_phylogenetic_tree/TerL.png"
    output:
        temp(touch(homedir + "/output/{sample}_4.TerL_phylogenetic_tree"))
    params:
        treedir = homedir + "/output/{sample}/TerL_tree"
    shell: '''# clean the temporary files in TerL tree
if [ -d {params.treedir} ];then rm -r {params.treedir};fi
'''
