# Snakemake script
clusterone_path = homedir + '/external_software/clusterone/cluster_one-1.0.jar'
runvc_script_dir = script_dir + '/run_vConTACT2'
rule gene2genome:
    input: homedir + '/output/{sample}/1.annotations/{sample}.faa'
    output: homedir + '/output/{sample}/5.vConTACT2_network/vcontact2/gene_to_genome.csv'
    log: log_dir + "/{sample}_run_vConTACT2.log"
    run:
        # 1 vConTACT gene2genome
        with open(output[0],'w') as f:
            for line in open(input[0]):
                if '>' in line:
                    name = line.split(' ',1)[0].replace('>', '')
                    f.write(f'{name},{wildcards.sample},none\n')

rule Diamond_blastp_raw:
    input:
        sample_faa = homedir + '/output/{sample}/1.annotations/{sample}.faa'
    output: 
        sample_dmnd = homedir + '/output/{sample}/5.vConTACT2_network/vcontact2/{sample}.dmnd',
        blp_f = homedir + '/output/{sample}/5.vConTACT2_network/vcontact2/{sample}_f.dm.tsv',
        blp_r = homedir + '/output/{sample}/5.vConTACT2_network/vcontact2/{sample}_r.dm.tsv',
        final_blp = homedir + '/output/{sample}/5.vConTACT2_network/vcontact2/{sample}_bidir.dm.tsv'
    params: 
        database_faa = f"{db_path}/{db_prefix}_vConTACT2_proteins.faa",
        vc2_db = f"{db_path}/{db_prefix}_vConTACT2_proteins.dmnd"
    threads: 60
    log: log_dir + "/{sample}_annotations.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 2 bidirectory blast to get vConTACT input blastp file
# the blastp cat front & reverse blastp outputs to final_blp
diamond blastp --query {input.sample_faa} --db {params.vc2_db} \\
    -o {output.blp_f} -p {threads} --sensitive --quiet >> {log}
diamond makedb --in {input.sample_faa} -p {threads} --db {output.sample_dmnd}
diamond blastp --query {params.database_faa} --db {output.sample_dmnd} \\
    -o {output.blp_r} -p {threads} --sensitive --quiet >> {log}
cat {output.blp_f} {output.blp_r} > {output.final_blp}
'''

rule vConTACT2:
    input: 
        g2g = homedir + '/output/{sample}/5.vConTACT2_network/vcontact2/gene_to_genome.csv',
        blp = homedir + '/output/{sample}/5.vConTACT2_network/vcontact2/{sample}_bidir.dm.tsv'
    output:
        network = homedir + '/output/{sample}/5.vConTACT2_network/c1.ntw',
        overview = homedir + '/output/{sample}/5.vConTACT2_network/genome_by_genome_overview.csv'
    params:
        wk_dir = homedir + '/output/{sample}/5.vConTACT2_network/vcontact2',
        vcontact2_db_blp = f"{db_path}/allVSall.dm.tsv",
        vcontact2_db = f"{db_path}/{db_prefix}_vConTACT2_proteins.dmnd",
        vcontact2_db_g2g = f"{db_path}/{db_prefix}_vConTACT2_gene_to_genome.csv"
    log: log_dir + "/{sample}_run_vConTACT2.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    threads: 60
    shell: '''# 3 vConTACT2 main protocol
# fmt input files, the sample head line was removed.
cat {params.vcontact2_db_blp} {input.blp} > {params.wk_dir}/vConTACT2_blastp.tsv
cat {params.vcontact2_db_g2g} {input.g2g} > {params.wk_dir}/vConTACT2_gene_to_genome.csv
# vcontact2 protocol
vcontact2 --blast-fp {params.wk_dir}/vConTACT2_blastp.tsv \\
    --proteins-fp {params.wk_dir}/vConTACT2_gene_to_genome.csv \\
    --db 'None' --pcs-mode MCL --vcs-mode ClusterONE -t {threads} \\
    --c1-bin {clusterone_path} --output-dir {params.wk_dir} >> {log}
# get final output
cp {params.wk_dir}/c1.ntw {output.network}
cp {params.wk_dir}/genome_by_genome_overview.csv {output.overview}
'''

rule vConTACT_visualize:
    input: 
        network = homedir + '/output/{sample}/5.vConTACT2_network/c1.ntw',
        overview = homedir + '/output/{sample}/5.vConTACT2_network/genome_by_genome_overview.csv',
        db_data = f'{db_path}/{db_prefix}_data.tsv'
    output: 
        small_network = homedir + '/output/{sample}/5.vConTACT2_network/small_c1.ntw',
        small_overview = homedir + '/output/{sample}/5.vConTACT2_network/small_genome_by_genome_overview.csv',
        small_database = homedir + '/output/{sample}/5.vConTACT2_network/small_data.tsv',
        graph = homedir + "/output/{sample}/5.vConTACT2_network/{sample}_vConTACT2.html"
    log: log_dir + "/{sample}_run_vConTACT2.log"
    threads: 60
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    params: outdir = homedir + '/output/{sample}/5.vConTACT2_network/graphanalyze_'
    shell: '''# 4 vConTACT visualization
python {runvc_script_dir}/smaller_vc_network.py -nwk {input.network} \\
    -ovv {input.overview} -db {input.db_data} -s {wildcards.sample} \\
    -onwk {output.small_network} -oovv {output.small_overview} \\
    -odb {output.small_database} >> {log}
if [ ! -d {params.outdir} ];then mkdir -p {params.outdir};fi
if [ -s {output.small_network} ];then
    python {runvc_script_dir}/graphanalyzer.py --graph {input.network} --csv {input.overview} \\
        --metas {input.db_data} --prefix {wildcards.sample} --suffix '' \\
        --view 2d -t {threads} -o {params.outdir} >> {log}
    cp {params.outdir}single-views_/{wildcards.sample}.html {output.graph}
    mv {params.outdir}single-views_ {params.outdir}
    mv {params.outdir}.pickle {params.outdir}
    mv {params.outdir}results_vcontact2_.xlsx {params.outdir}
    mv {params.outdir}results_vcontact2_.csv {params.outdir}
    mv {params.outdir}graphanalyzer_.log {params.outdir}
    mv {params.outdir}custom_taxonomy_table__mqc.txt {params.outdir}
    mv {params.outdir}csv_edit_.xlsx {params.outdir}
else
    echo "The phage is singleton in INPHARED database."
    touch {output.graph}
fi
'''

rule vcontact_all:
    input:
        network = homedir + '/output/{sample}/5.vConTACT2_network/c1.ntw',
        overview = homedir + '/output/{sample}/5.vConTACT2_network/genome_by_genome_overview.csv',
        small_network = homedir + '/output/{sample}/5.vConTACT2_network/small_c1.ntw',
        small_overview = homedir + '/output/{sample}/5.vConTACT2_network/small_genome_by_genome_overview.csv',
        small_database = homedir + '/output/{sample}/5.vConTACT2_network/small_data.tsv',
        graph = homedir + "/output/{sample}/5.vConTACT2_network/{sample}_vConTACT2.html",
        db_data = f'{db_path}/{db_prefix}_data.tsv'
    output:
        touch(temp(homedir + '/output/{sample}_5.vConTACT2_network')),
        db_data = homedir + '/output/{sample}/5.vConTACT2_network/' + f'{db_prefix}_data.tsv'
    params:
        workdir = homedir + "/output/{sample}/5.vConTACT2_network/vcontact2"
    shell: '''# 5 cleanup the tempory
if [ -d {params.workdir} ];then rm -r {params.workdir};fi
cp {input.db_data} {output.db_data}
'''
