# 2.4.1 vConTACT gene2genome
rule gene2genome:
    input: 'output/{sample}/{sample}.faa'
    output: 'output/{sample}/vcontact2/gene_to_genome.csv'
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
        vcontact2_db_blp = f"{db_path}allVSall.dm.tsv",
        vcontact2_db = f"{db_path}{db_prefix}_vConTACT2_proteins.dmnd",
        vcontact2_db_g2g = f"{db_path}{db_prefix}_vConTACT2_gene_to_genome.csv"
    conda: "envs/phagesnake.yaml"
    threads: 60
    shell: '''cat {params.vcontact2_db_blp} {input.blp} > \\
    {params.wk_dir}/vConTACT2_blastp.tsv
sed -i "1d" {input.g2g} # delete table header
cat {params.vcontact2_db_g2g} {input.g2g} > \\
    {params.wk_dir}/vConTACT2_gene_to_genome.csv

vcontact2 --blast-fp {params.wk_dir}/vConTACT2_blastp.tsv \\
    --proteins-fp {params.wk_dir}/vConTACT2_gene_to_genome.csv \\
    --db 'None' --rel-mode Diamond --pcs-mode MCL \\
    --vcs-mode ClusterONE --output-dir {params.wk_dir} -t {threads}

cp {params.wk_dir}/c1.ntw {output.network}
cp {params.wk_dir}/genome_by_genome_overview.csv {output.overview}
'''


# 2.4.3 vConTACT visualization
# rule vConTACT_visualize:
#     input: 
#         network = 'output/{sample}/vcontact2/c1.ntw',
#         overview = 'output/{sample}/vcontact2/genome_by_genome_overview.csv'
#     output: "output/{sample}/{sample}_vConTACT.pdf"
#     conda: "envs/phagesnake.yaml"
#     shell: ''''''

