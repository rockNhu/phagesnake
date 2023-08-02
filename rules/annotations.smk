# Snakemake script
# 2.1 prodigal annotation
rule Prodigal:
    input: fmt_fna_dir +"/{sample}.fasta"
    output: 
        faa = "output/{sample}/{sample}.faa",
        gff = "output/{sample}/{sample}.gff"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    log: f"{log_dir}/" + "{sample}_annotations.log"
    shell: "prodigal -i {input} -a {output.faa} -f gff -o {output.gff} -c -q -p meta > {log}"

# 2.1.1 EggNOG annotation
rule EggNOG:
    input: "output/{sample}/{sample}.faa"
    output: "output/{sample}/eggnog/{sample}.emapper.annotations"
    threads: 50
    params:
        out_dir = 'output/{sample}/eggnog'
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''if [ ! -d {params.out_dir} ];then mkdir -p {params.out_dir};fi
emapper.py -i {input} -o {wildcards.sample} --cpu {threads} \\
    --output_dir {params.out_dir} --override --go_evidence non-electronic >> {log}
'''

# 2.2.1 Diamond blastp All protein to find the annotation
rule Diamond_blastp:
    input: 'output/{sample}/{sample}.faa'
    output: 'output/{sample}/{sample}.dm.tsv'
    params: 
        vc2_db = f"{db_path}/{db_prefix}_vConTACT2_proteins.dmnd"
    threads: 60
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''diamond blastp --query {input} --db {params.vc2_db} -o {output} \\
    -p {threads} --quiet --sensitive --outfmt 6 qseqid sseqid pident qcovhsp >> {log}
'''

# 2.2.2 parse the blastp All output
rule filter_All_BLASTp:
    input: 'output/{sample}/{sample}.dm.tsv'
    output: 'output/{sample}/blastp_fmt.tsv'
    params: 
        totalname_dict = f'{db_path}/{db_prefix}_vConTACT2_proteins_totalname.pydict'
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
        fa = fmt_fna_dir + "/{sample}.fasta",
        faa = "output/{sample}/{sample}.faa",
        blastp = "output/{sample}/blastp_fmt.tsv",
        eggnog = "output/{sample}/eggnog/{sample}.emapper.annotations"
    output: "output/{sample}/{sample}.gbk"
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''python {script_dir}/make_final_gbk.py \\
    -f {input.fa} -a {input.faa} \\
    -b {input.blastp} -e {input.eggnog} \\
    -o {output} >> {log}
'''

# 2.3.1 genome_viewer
rule genome_visualize:
    input: "output/{sample}/{sample}.gbk"
    output: 
        png_out = "output/{sample}/{sample}.png",
        svg_out = "output/{sample}/{sample}.svg"
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''python {script_dir}/plot_arrow.py -i {input} -o {output.png_out} >> {log}
python {script_dir}/plot_arrow.py -i {input} -o {output.svg_out} >> {log}
'''

# abricate protocol
rule abr:
    input: fmt_fna_dir +"/{sample}.fasta"
    output: 
        card_out = "output/{sample}/ABRicate/CARD_out.tab",
        vfdb_out = "output/{sample}/ABRicate/VFDB_out.tab",
        resf_out = "output/{sample}/ABRicate/ResF_out.tab",
        arga_out = "output/{sample}/ABRicate/ArgA_out.tab",
        ecoh_out = "output/{sample}/ABRicate/EcoH_out.tab",
        ecov_out = "output/{sample}/ABRicate/Ecoli_VF_out.tab",
        megares_out = "output/{sample}/ABRicate/MegaRes_out.tab",
        ncbi_out = "output/{sample}/ABRicate/NCBI_out.tab",
        abr_check = "output/{sample}/abr_check.tsv"
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    threads: 8
    params:
        out_dir = "output/{sample}/ABRicate"
    shell: '''if [ ! -d {params.out_dir} ];then mkdir -p {params.out_dir};fi
abricate --quiet --nopath --db card {input} 1> {output.card_out} 2>> {log}
abricate --quiet --nopath --db vfdb {input} 1> {output.vfdb_out} 2>> {log}
abricate --quiet --nopath --db resfinder {input} 1> {output.resf_out} 2>> {log}
abricate --quiet --nopath --db argannot {input} 1> {output.arga_out} 2>> {log}
abricate --quiet --nopath --db ecoh {input} 1> {output.ecoh_out} 2>> {log}
abricate --quiet --nopath --db ecoli_vf {input} 1> {output.ecov_out} 2>> {log}
abricate --quiet --nopath --db megares {input} 1> {output.megares_out} 2>> {log}
abricate --quiet --nopath --db ncbi {input} 1> {output.ncbi_out} 2>> {log}
python {script_dir}/count_final.py -i {params.out_dir} -o {output.abr_check} >> {log}
'''