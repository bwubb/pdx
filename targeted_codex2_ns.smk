import errno
import os

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno==errno.EEXIST and os.path.isdir(path):
            pass

with open(config['project']['sample_list'],'r') as i:
    SAMPLES=i.read().splitlines()

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

for sample in SAMPLES:
    mkdir_p('logs/cluster/%s' % sample)

rule all:
    input:
        f"data/work/{config['project']['name']}/codex2/segments.txt"

rule CODEX2_ns:
    input:
        [BAMS[sample] for sample in SAMPLES]
    output:
        "data/work/{project}/codex2/segments.txt"
    params:
        project="data/work/{project}/codex2",
        bed=config['resources']['targets_bed'],
        normals=config['project']['normal_list']
    script:
        "codex2_ns_targeted.R"
