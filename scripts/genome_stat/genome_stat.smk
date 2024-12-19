# Snakemake script
# The genomic statistic sub-workflow
stat_script_dir = script_dir + '/genome_stat'

rule get_orf_list:
    input: 
        faas = expand(homedir + "/output/{sample}/1.annotations/{sample}.faa",sample=Samples)
    output:
        orf = temp(homedir + "/output/orf.list")
    log: f"{log_dir}/genome_stat.log"
    run:
        # 3.1 get orf.list
        import os

        out = open(output.orf,'w')
        for faa in input.faas:
            name = os.path.basename(faa).rsplit(".")[0]
            count = sum(1 for line in open(faa) if '>' in line)
            out.write(f'{name}:{count}\n')

rule statistics_genome:
    input:
        orf = homedir + "/output/orf.list",
        fna_dir = fmt_fna_dir
    output: f"{homedir}/output/2.seq_info{start_time}.tsv"
    log: f"{log_dir}/genome_stat.log"
    #conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''# 3.2 statistic genome infomation
python {stat_script_dir}/get_fasta_info.py \\
    -i {input.fna_dir} -orf {input.orf} -o {output} >> {log}
'''
