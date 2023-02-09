# 3.1 get orf.list
rule get_orf_list:
    input: 
        faas = expand("output/{sample}/{sample}.faa",sample=Samples)
    output: "orf.list"
    run:
        import os

        out = open(output[0],'w')
        for faa in input.faas:
            name = os.path.basename(faa).rsplit(".")[0]
            count = sum(1 for line in open(faa) if '>' in line)
            out.write(f'{name}:{count}\n')


# 3.2 statistic genome infomation
rule statistics_genome:
    input:
        fna_dir = fna_dir,
        orf = "orf.list"
    output: "seq_info.tsv"
    params: 
        scripts = script_dir
    conda: "envs/phagesnake.yaml"
    shell: '''python {params.scripts}/get_fasta_info.py \\
    -i {input.fna_dir} -orf {input.orf} -o .
'''