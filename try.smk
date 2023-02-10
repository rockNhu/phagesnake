configfile: '../config.yaml'
db_path = config['db_path']
db_prefix = config['db_prefix']

rule echo_x:
    output: temp(touch("ok"))
    shell:''' echo {db_path}
echo {db_prefix}
'''