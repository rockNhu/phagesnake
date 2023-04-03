# Snakemake script
# 3.1 get orf.list
rule get_orf_list:
    input: 
        faas = expand("output/{sample}/{sample}.faa",sample=Samples)
    output:
        orf = "orf.list"
    log: f"{log_dir}/genome_stat.log"
    run:
        import os

        out = open(output.orf,'w')
        for faa in input.faas:
            name = os.path.basename(faa).rsplit(".")[0]
            count = sum(1 for line in open(faa) if '>' in line)
            out.write(f'{name}:{count}\n')


# 3.2 statistic genome infomation
rule statistics_genome:
    input:
        orf = "orf.list"
    output: protected("seq_info.tsv")
    log: f"{log_dir}/genome_stat.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    shell: '''python {script_dir}/get_fasta_info2.py \\
    -i {fmt_fna_dir} -orf {input.orf} -o . >> {log}
'''