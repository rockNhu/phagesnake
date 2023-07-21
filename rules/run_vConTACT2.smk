# Snakemake script

# 2.4.1 vConTACT gene2genome
rule gene2genome:
    input: 'output/{sample}/{sample}.faa'
    output: 'output/{sample}/vcontact2/gene_to_genome.csv'
    log: f"{log_dir}/" + "{sample}_run_vConTACT2.log"
    run:
        with open(output[0],'w') as f:
            for line in open(input[0]):
                if '>' in line:
                    name = line.split(' ',1)[0].replace('>', '')
                    f.write(f'{name},{wildcards.sample},none\n')

#TODO: 2.2.1 the fmt checking
rule Diamond_blastp_raw:
    input:
        sample_faa = 'output/{sample}/{sample}.faa'
    output: 
        sample_dmnd = 'output/{sample}/vcontact2/{sample}.dmnd',
        blp_f = 'output/{sample}/vcontact2/{sample}_f.dm.tsv',
        blp_r = 'output/{sample}/vcontact2/{sample}_r.dm.tsv',
        final_blp = 'output/{sample}/vcontact2/{sample}_bidir.dm.tsv'
    params: 
        database_faa = f"{db_path}/{db_prefix}_vConTACT2_proteins.faa",
        vc2_db = f"{db_path}/{db_prefix}_vConTACT2_proteins.dmnd"
    threads: 60
    log: f"{log_dir}/" + "{sample}_annotations.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# the blastp cat front & reverse blastp outputs to final_blp
diamond blastp --query {input.sample_faa} --db {params.vc2_db} \\
    -o {output.blp_f} -p {threads} --sensitive --quiet >> {log}
diamond makedb --in {input.sample_faa} -p {threads} --db {output.sample_dmnd}
diamond blastp --query {params.database_faa} --db {output.sample_dmnd} \\
    -o {output.blp_r} -p {threads} --sensitive --quiet >> {log}
cat {output.blp_f} {output.blp_r} > {output.final_blp}
'''

# 2.4.2 vConTACT2
rule vConTACT2:
    input: 
        g2g = 'output/{sample}/vcontact2/gene_to_genome.csv',
        blp = 'output/{sample}/vcontact2/{sample}_bidir.dm.tsv'
    output:
        network = 'output/{sample}/c1.ntw',
        overview = 'output/{sample}/genome_by_genome_overview.csv'
    params:
        wk_dir = 'output/{sample}/vcontact2',
        vcontact2_db_blp = f"{db_path}/allVSall.dm.tsv",
        vcontact2_db = f"{db_path}/{db_prefix}_vConTACT2_proteins.dmnd",
        vcontact2_db_g2g = f"{db_path}/{db_prefix}_vConTACT2_gene_to_genome.csv"
    log: f"{log_dir}/" + "{sample}_run_vConTACT2.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    threads: 60
    shell: '''# fmt input files, the sample head line was removed.
cat {params.vcontact2_db_blp} {input.blp} > {params.wk_dir}/vConTACT2_blastp.tsv
cat {params.vcontact2_db_g2g} {input.g2g} > {params.wk_dir}/vConTACT2_gene_to_genome.csv

vcontact2 --blast-fp {params.wk_dir}/vConTACT2_blastp.tsv \\
    --proteins-fp {params.wk_dir}/vConTACT2_gene_to_genome.csv \\
    --db 'None' --pcs-mode MCL --vcs-mode ClusterONE -t {threads} \\
    --output-dir {params.wk_dir} >> {log}

cp {params.wk_dir}/c1.ntw {output.network}
cp {params.wk_dir}/genome_by_genome_overview.csv {output.overview}
'''

# 2.4.3 vConTACT visualization
rule vConTACT_visualize_test:
    input: 
        network = 'output/{sample}/c1.ntw',
        overview = 'output/{sample}/genome_by_genome_overview.csv',
        db_data = f'{db_path}/{db_prefix}_data.tsv'
    output: 
        small_network = 'output/{sample}/small_c1.ntw',
        small_overview = 'output/{sample}/small_genome_by_genome_overview.csv',
        small_database = 'output/{sample}/small_data.tsv',
        graph = "output/{sample}/{sample}.html"
    log: f"{log_dir}/" + "{sample}_run_vConTACT2.log"
    threads: 60
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    params: outdir = 'output/{sample}/graphanalyze'
    shell: '''
python {script_dir}/smaller_vc_network.py -nwk {input.network} \\
    -ovv {input.overview} -db {input.db_data} -s {wildcards.sample} \\
    -onwk {output.small_network} -oovv {output.small_overview} \\
    -odb {output.small_database} >> {log}
if [ ! -d {params.outdir} ];mkdir -p {params.outdir};fi
if [ -s {output.small_network} ];then
    python {script_dir}/graphanalyzer.py --graph {input.network} --csv {input.overview} \\
        --metas {input.db_data} --prefix {wildcards.sample} --suffix '' \\
        --view 2d -t {threads} -o {params.outdir} >> {log}
    cp {params.outdir}/*.html {output.graph}
else
    echo "The phage is singleton in INPHARED database."
    touch {output.graph}
fi
'''