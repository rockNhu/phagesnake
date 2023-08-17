# Snakemake script
# The annotations sub-workflow

rule Prodigal:
    input: fmt_fna_dir +"/{sample}.fasta"
    output: 
        faa = "output/{sample}/{sample}.faa",
        gff = "output/{sample}/{sample}.gff"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    log: f"{log_dir}/" + "{sample}_annotations.log"
    shell: '''# 2.1 prodigal annotation
prodigal -i {input} -a {output.faa} -f gff -o {output.gff} -c -q -p meta > {log}
'''

rule EggNOG:
    input: "output/{sample}/{sample}.faa"
    output: "output/{sample}/eggnog/{sample}.emapper.annotations"
    threads: 50
    params:
        out_dir = 'output/{sample}/eggnog'
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2.1.1 EggNOG annotation
if [ ! -d {params.out_dir} ];then mkdir -p {params.out_dir};fi
emapper.py -i {input} -o {wildcards.sample} --cpu {threads} \
--output_dir {params.out_dir} --override --go_evidence non-electronic >> {log}
'''

rule Diamond_blastp:
    input: 'output/{sample}/{sample}.faa'
    output: temp('output/{sample}/{sample}.dm.tsv')
    params: 
        vc2_db = f"{db_path}/{db_prefix}_vConTACT2_proteins.dmnd"
    threads: 60
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2.2.1 Diamond blastp All protein to find the annotation
diamond blastp --query {input} --db {params.vc2_db} -o {output} \
-p {threads} --quiet --sensitive --outfmt 6 qseqid sseqid pident qcovhsp >> {log}
'''

rule filter_All_BLASTp:
    input: tsv = 'output/{sample}/{sample}.dm.tsv'
    output: fmt_blastp = 'output/{sample}/blastp_fmt.tsv'
    params: 
        nucl_totalname_dict = f'{db_path}/{db_prefix}_genomes_totalname.pydict',
        prot_totalname_dict = f'{db_path}/{db_prefix}_vConTACT2_proteins_totalname.pydict',
        idt = 70,
        cov = 70
    log: f"{log_dir}/" + "{sample}_annotations.log"
    shell: '''# 2.2.2 parse the blastp All output
python {script_dir}/get_neibour_prot.py \
-i {input.tsv} -o_fmt {output.fmt_blastp} \
-nucl_id_db {params.nucl_totalname_dict} -prot_id_db {params.prot_totalname_dict} \
-idt {params.idt} -cov {params.cov} >> {log}
'''

rule final_gbk:
    input:
        fa = fmt_fna_dir + "/{sample}.fasta",
        faa = "output/{sample}/{sample}.faa",
        blastp = "output/{sample}/blastp_fmt.tsv",
        eggnog = "output/{sample}/eggnog/{sample}.emapper.annotations"
    output: "output/{sample}/{sample}.gbk"
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2.3 fasta + Prodigal + EggNOG + blastp -> gbk
python {script_dir}/make_final_gbk.py \
-f {input.fa} -a {input.faa} \
-b {input.blastp} -e {input.eggnog} \
-o {output} >> {log}
'''

rule genome_visualize:
    input: "output/{sample}/{sample}.gbk"
    output: 
        png_out = "output/{sample}/{sample}.png",
        svg_out = "output/{sample}/{sample}.svg"
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2.3.1 genome_viewer
python {script_dir}/plot_arrow.py -i {input} -o {output.png_out} >> {log}
python {script_dir}/plot_arrow.py -i {input} -o {output.svg_out} >> {log}
'''

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
    shell: '''# 2.4 abricate annotation
if [ ! -d {params.out_dir} ];then mkdir -p {params.out_dir};fi
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