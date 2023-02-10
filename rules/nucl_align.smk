# Snakemake script
# 1.1 blastn to find the similar species genomes
rule MMseqs_blastn:
    input:
        fna_dir + "/{sample}.fasta"
    output:
        'output/{sample}/blastn.tsv'
    threads: 20
    conda: "envs/phagesnake.yaml"
    params:
        db_path = f'{db_path}{db_prefix}_genomes.fa'
    shell: '''mmseqs easy-search --search-type 3 {input} \\
    {params.db_path} {output} output/{wildcards.sample}/tmp \\
    --threads {threads} --format-output "query,target,pident,qcov"
# clean up the temporary
rm -rf output/{wildcards.sample}/tmp
'''


# 1.2 parse the blastn output
rule catch_nucl_neibours:
    input: 
        tsv = 'output/{sample}/blastn.tsv'
    output: 
        l = 'output/{sample}/blastn.list'
    run: 
        to_extract = dict()
        for line in open(input.tsv):
            content = line.strip('\n').split('\t') # input type is tsv
            acc = content[1] # the second column is accession number
            idt = float(content[2]) # the third column is fraction identity
            cov = float(content[3]) # the fourth column is fraction coverage
            if acc not in to_extract:
                to_extract[acc] = [[idt],[cov]]
            else:
                to_extract[acc][0].append(idt) # the idtentity list
                to_extract[acc][1].append(cov) # the coverage list
        # remove idt and cov away from one genus
        final_acc = []
        for acc, values in to_extract.items():
            final_idt = sum(values[0]) / len(values[0]) # the idt is avg
            final_cov = sum(values[1]) * 100 # sum up cov and convert it to perc
            if final_idt >= 75 and final_cov >= 75:
                final_acc.append(acc)
        # output to name list
        with open(output.l,'w') as f:
            f.writelines(i  + '\n' for i in final_acc)


# 1.3 catch the neibour fna and seperate neibours wtih format to "n_output"
rule catch_neibours_fna:
    input: 
        fa = fna_dir + "/{sample}.fasta",
        ex_list = 'output/{sample}/blastn.list'
    output: 
        directory("output/{sample}/n_output")
    conda: "envs/phagesnake.yaml"
    params:
        script_dir = script_dir,
        totalname_dict = f'{db_path}{db_prefix}_genomes_totalname.pydict',
        nameseq_dict = f'{db_path}{db_prefix}_genomes_nameseq.pydict'
    shell: '''python {params.script_dir}/get_seqs_from_dict.py -i {input.ex_list} \\
    -o {output} -ns {params.nameseq_dict} -tn {params.totalname_dict} \\
    --seperate
# add self fna
cp {input.fa} {output}
'''


# 1.4 do ANI using the pyani, the method is ANIb
rule pyANI:
    input: 
        to_check = 'output/{sample}/blastn.list',
        nd = "output/{sample}/n_output" # blastn out dir
    output:
        directory("output/{sample}/ANI_output")
    conda: "envs/phagesnake.yaml"
    shell: '''if [ -s {input.to_check} ];then
    average_nucleotide_identity.py -i {input.nd} -o {output} -m ANIb -g
else
    echo "Error: empty neibours with {wildcards.sample}"
    mkdir -p {output}
fi
'''
