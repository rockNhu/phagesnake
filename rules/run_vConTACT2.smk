# Snakemake script
# 2.4.1 vConTACT gene2genome
rule gene2genome:
    input: 'output/{sample}/{sample}.faa'
    output: 'output/{sample}/vcontact2/gene_to_genome.csv'
    log: f"{log_dir}/" + "{sample}_run_vConTACT2.log"
    run:
        with open(output[0],'w') as f:
            f.write('protein_id,contig_id,keywords\n')
            for line in open(input[0]):
                if '>' in line:
                    name = line.split(' ',1)[0].replace('>', '')
                    f.write(f'{name},{wildcards.sample},None_provided\n')


# 2.4.2 vConTACT2
rule vConTACT2:
    input: 
        faa = 'output/{sample}/{sample}.faa',
        g2g = 'output/{sample}/vcontact2/gene_to_genome.csv',
        blp = 'output/{sample}/{sample}.dm.tsv'
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
    shell: '''cat {params.vcontact2_db_blp} {input.blp} > \\
    {params.wk_dir}/vConTACT2_blastp.tsv
sed -i "1d" {input.g2g} # delete table header
cat {params.vcontact2_db_g2g} {input.g2g} > \\
    {params.wk_dir}/vConTACT2_gene_to_genome.csv

vcontact2 --blast-fp {params.wk_dir}/vConTACT2_blastp.tsv \\
    --proteins-fp {params.wk_dir}/vConTACT2_gene_to_genome.csv \\
    --db 'None' --rel-mode Diamond --pcs-mode MCL \\
    --vcs-mode ClusterONE --output-dir {params.wk_dir} -t {threads} >> {log}

cp {params.wk_dir}/c1.ntw {output.network}
cp {params.wk_dir}/genome_by_genome_overview.csv {output.overview}
'''


# 2.4.3 vConTACT visualization
rule vConTACT_visualize:
    input: 
        network = 'output/{sample}/c1.ntw',
        overview = 'output/{sample}/genome_by_genome_overview.csv',
        db_data = f'{db_path}/{db_prefix}_data.tsv'
    output: 
        small_network = temp('output/{sample}/small_c1.ntw'),
        small_overview = temp('output/{sample}/small_genome_by_genome_overview.csv'),
        small_database = temp('output/{sample}/small_data.tsv'),
        big_graph = directory("output/{sample}/{sample}_vConTACT_plot"),
        small_graph = directory("output/{sample}/{sample}_small_vConTACT_plot")
    log: f"{log_dir}/" + "{sample}_run_vConTACT2.log"
    threads: 60
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''
python {script_dir}/smaller_vc_network.py -nwk {input.network} \\
    -ovv {input.overview} -db {input.db_data} -s {wildcards.sample} \\
    -onwk {output.small_network} -oovv {output.small_overview} \\
    -odb {output.small_database} >> {log}
if [ -s {output.small_network} ];then
    # the big graph
    python {script_dir}/graphanalyzer.py --graph {input.network} --csv {input.overview} \\
        --metas {input.db_data} --prefix {wildcards.sample} --view 2d -t {threads} \\ 
        -o {output.big_graph} >> {log}
    # the small graph
    python {script_dir}/graphanalyzer.py --graph {input.network} --csv {input.overview} \\
        --metas {input.db_data} --prefix {wildcards.sample} --view 2d -t {threads} \\
        -o {output.small_graph} >> {log}
else
    echo "The phage is singleton in INPHARED database."
    mkdir -p {output.big_graph}
    mkdir -p {output.small_graph}
'''

