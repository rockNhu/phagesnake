# Snakemake script
# The annotations sub-workflow

anno_script_dir = script_dir + '/annotations'

rule Prodigal:
    input: fmt_fna_dir + "/{sample}.fasta"
    output: 
        faa = homedir + "/output/{sample}/2.annotations/{sample}.faa",
        gff = homedir + "/output/{sample}/2.annotations/{sample}.gff"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    log: log_dir + "/{sample}_annotations.log"
    shell: '''# 2.1 prodigal annotation
prodigal -i {input} -a {output.faa} -f gff -o {output.gff} -c -q -p meta > {log}
'''

rule EggNOG:
    input: homedir + "/output/{sample}/2.annotations/{sample}.faa"
    output: homedir + "/output/{sample}/2.annotations/eggnog/{sample}.emapper.annotations"
    threads: 50
    params:
        out_dir = homedir + '/output/{sample}/2.annotations/eggnog'
    log: log_dir + "/{sample}_annotations.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2.1.1 EggNOG annotation
if [ ! -d {params.out_dir} ];then mkdir -p {params.out_dir};fi
emapper.py --data_dir {db_path} -i {input} -o {wildcards.sample} --cpu {threads} \
--output_dir {params.out_dir} --override --go_evidence non-electronic >> {log}
'''

rule Diamond_blastp:
    input: homedir + '/output/{sample}/2.annotations/{sample}.faa'
    output: temp(homedir + '/output/{sample}/2.annotations/{sample}.dm.tsv')
    params: 
        vc2_db = f"{db_path}/{db_prefix}_vConTACT2_proteins.dmnd"
    threads: 60
    log: log_dir + "/{sample}_annotations.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2.2.1 Diamond blastp All protein to find the annotation
diamond blastp --query {input} --db {params.vc2_db} -o {output} \
-p {threads} --quiet --sensitive --outfmt 6 qseqid sseqid pident qcovhsp >> {log}
'''

rule filter_All_BLASTp:
    input: tsv = homedir + '/output/{sample}/2.annotations/{sample}.dm.tsv'
    output: fmt_blastp = homedir + '/output/{sample}/2.annotations/blastp_fmt.tsv'
    params: 
        nucl_totalname_dict = f'{db_path}/{db_prefix}_genomes_totalname.pydict',
        prot_totalname_dict = f'{db_path}/{db_prefix}_vConTACT2_proteins_totalname.pydict',
        idt = 70,
        cov = 70
    log: log_dir + "/{sample}_annotations.log"
    shell: '''# 2.2.2 parse the blastp All output
python {anno_script_dir}/get_neibour_prot.py \
-i {input.tsv} -o_fmt {output.fmt_blastp} \
-nucl_id_db {params.nucl_totalname_dict} -prot_id_db {params.prot_totalname_dict} \
-idt {params.idt} -cov {params.cov} >> {log}
'''

rule final_gbk:
    input:
        fa = fmt_fna_dir + "/{sample}.fasta",
        faa = homedir + "/output/{sample}/2.annotations/{sample}.faa",
        blastp = homedir + "/output/{sample}/2.annotations/blastp_fmt.tsv",
        eggnog = homedir + "/output/{sample}/2.annotations/eggnog/{sample}.emapper.annotations"
    output: homedir + "/output/{sample}/2.annotations/{sample}.gbk"
    log: log_dir + "/{sample}_annotations.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2.3 fasta + Prodigal + EggNOG + blastp -> gbk
python {anno_script_dir}/make_final_gbk.py \
-f {input.fa} -a {input.faa} \
-b {input.blastp} -e {input.eggnog} \
-o {output} >> {log}
'''

rule genome_visualize:
    input: homedir + "/output/{sample}/2.annotations/{sample}.gbk"
    output: 
        png_out = homedir + "/output/{sample}/2.annotations/{sample}.png",
        svg_out = homedir + "/output/{sample}/2.annotations/{sample}.svg"
    log: log_dir + "/{sample}_annotations.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2.3.1 genome_viewer
python {anno_script_dir}/plot_arrow.py -i {input} -o {output.png_out} >> {log}
python {anno_script_dir}/plot_arrow.py -i {input} -o {output.svg_out} >> {log}
'''

rule abr:
    input: fmt_fna_dir +"/{sample}.fasta"
    output: 
        abr_check = homedir + "/output/{sample}/3.secuity_check/abr_check.tsv",
        abr_dir = directory(homedir + "/output/{sample}/3.secuity_check/ABRicate")
    log: log_dir + "/{sample}_abr_check.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    threads: 8
    shell: '''# 2.4 abricate annotation
if [ ! -d {output.abr_dir} ];then mkdir -p {output.abr_dir};fi
abricate --quiet --nopath --db card {input} 1> {output.abr_dir}/CARD_out.tab 2>> {log}
abricate --quiet --nopath --db vfdb {input} 1> {output.abr_dir}/VFDB_out.tab 2>> {log}
abricate --quiet --nopath --db resfinder {input} 1> {output.abr_dir}/ResF_out.tab 2>> {log}
abricate --quiet --nopath --db argannot {input} 1> {output.abr_dir}/ArgA_out.tab 2>> {log}
abricate --quiet --nopath --db ecoh {input} 1> {output.abr_dir}/EcoH_out.tab 2>> {log}
abricate --quiet --nopath --db ecoli_vf {input} 1> {output.abr_dir}/Ecoli_VF_out.tab 2>> {log}
#abricate --quiet --nopath --db megares {input} 1> {output.abr_dir}/MegaRes_out.tab 2>> {log}
abricate --quiet --nopath --db ncbi {input} 1> {output.abr_dir}/NCBI_out.tab 2>> {log}
python {anno_script_dir}/count_final.py -i {output.abr_dir} -o {output.abr_check} >> {log}
'''

rule annotation_all:
    input:
        png_out = homedir + "/output/{sample}/2.annotations/{sample}.png",
        svg_out = homedir + "/output/{sample}/2.annotations/{sample}.svg",
        gbk = homedir + "/output/{sample}/2.annotations/{sample}.gbk",
        gff = homedir + "/output/{sample}/2.annotations/{sample}.gff",
        faa = homedir + "/output/{sample}/2.annotations/{sample}.faa",
        blastp = homedir + "/output/{sample}/2.annotations/blastp_fmt.tsv",
        eggnog_out = homedir + "/output/{sample}/2.annotations/eggnog/{sample}.emapper.annotations",
        abr_dir = homedir + "/output/{sample}/3.secuity_check/ABRicate",
        abr_check = homedir + "/output/{sample}/3.secuity_check/abr_check.tsv"
    output:
        temp(touch(homedir + "/output/{sample}_2.annotations-clean")),
        temp(touch(homedir + "/output/{sample}_3.secuity_check-clean"))
