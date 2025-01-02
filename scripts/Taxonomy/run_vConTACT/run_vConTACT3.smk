# Snakemake script
runvc_script_dir = script_dir + '/run_vConTACT'

rule vConTACT3:
    input: 
        fa =  fmt_fna_dir + '/{sample}.fasta'
    output:
        homedir + '/output/{sample}/5.vConTACT_network/Done-vc3'
    params:
        wk_dir = homedir + '/output/{sample}/5.vConTACT_network/vcontact3',
        db_path = db_path
    log: log_dir + "/{sample}_run_vConTACT3.log"
    conda: f"{Conda_env_dir}/phagesnake.yaml"
    threads: 60
    shell: '''# 3 vConTACT2 main protocol
vcontact3 run --nucleotide {input.fa} --db-path {params.db_path} \\
    --db_domain "prokaryotes" -t {threads} --output {params.wk_dir} >> {log}
touch {output}
'''

rule vcontact_all:
    input:
        homedir + '/output/{sample}/5.vConTACT_network/Done-vc3'
    output:
        touch(temp(homedir + '/output/{sample}_5.vConTACT_network'))
