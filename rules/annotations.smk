# Snakemake script
# 2.1 prodigal annotation
rule Prodigal:
    input: "fna_files/{sample}.fasta"
    output: "output/{sample}/{sample}.faa"
    conda: "envs/phagesnake.yaml"
    shell: "prodigal -i {input} -a {output} -c -q"


# 2.1.1 EggNOG annotation
rule EggNOG:
    input: "output/{sample}/{sample}.faa"
    output: "output/{sample}/eggnog/{sample}.emapper.annotations"
    threads: 50
    params:
        out_dir = 'output/{sample}/eggnog'
    conda: "envs/phagesnake.yaml"
    shell: '''mkdir -p {params.out_dir}
emapper.py -i {input} -o {wildcards.sample} --cpu {threads} \\
    --output_dir {params.out_dir} --go_evidence non-electronic 
'''


# 2.2.1 Diamond blastp All protein to find the annotation
rule Diamond_blastp:
    input: 'output/{sample}/{sample}.faa'
    output: 'output/{sample}/{sample}.dm.tsv'
    params: 
        vc2_db = f"{db_path}{db_prefix}_vConTACT2_proteins.dmnd"
    conda: "envs/phagesnake.yaml"
    shell: '''#
diamond blastp --query {input} --db {params.vc2_db} -o {output} \\
-p 60 --sensitive --outfmt 6 qseqid sseqid pident qcovhsp
'''


# 2.2.2 parse the blastp All output
rule filter_All_BLASTp:
    input: 'output/{sample}/{sample}.dm.tsv'
    output: 'output/{sample}/blastp_fmt.tsv'
    params: 
        totalname_dict = f'{db_path}{db_prefix}_vConTACT2_proteins_totalname.pydict'
    run:
        def is_unknown(product):
            unknow_words = [
                'hypothetical','Unknow','unknow',
                'Uncharactrized','uncharactrized','gp','Gp'
            ]
            return any(un in product for un in unknow_words)


        # load the id:totalname
        with open(params.totalname_dict) as f:
            total_name = eval(f.read())

        out_line = ''
        for line in open(input[0]):
            contents = line.strip('\n').split('\t')
            name = contents[0]
            target = contents[1]
            product = total_name[target].split(' ',1)[1]
            idt = float(contents[2])
            qcov = float(contents[3])
            # the filter
            if not is_unknown(product) or qcov < 70 or idt < 70:
                out_line += f'{name}\t{product}\t{qcov}\t{idt}\t{target}\n'

        with open(output[0],'w') as f:
            f.write(out_line)


# 2.3 fasta + Prodigal + EggNOG + blastp -> gbk
rule final_gbk:
    input:
        fa = fna_dir + "/{sample}.fasta",
        faa = "output/{sample}/{sample}.faa",
        blastp = "output/{sample}/blastp_fmt.tsv",
        eggnog = "output/{sample}/eggnog/{sample}.emapper.annotations"
    output: "output/{sample}/{sample}.gbk"
    conda: "envs/phagesnake.yaml"
    shell: '''python {script_dir}/make_final_gbk.py \\
    -f {input.fa} -a {input.faa} \\
    -b {input.blastp} -e {input.eggnog} \\
    -o {output}
'''


# 2.3.1 genome_viewer
rule genome_visualize:
    input: "output/{sample}/{sample}.gbk"
    output: "output/{sample}/{sample}.png"
    conda: "envs/phagesnake.yaml"
    shell: "python {script_dir}/plot_arrow.py -i {input} -o {output}"
