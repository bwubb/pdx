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

with open(config['project']['bam_list'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

for sample in SAMPLES:
    mkdir_p('logs/cluster/%s' % sample)

rule all:
    input:
        f"data/work/{config['project']['name']}/codex2/coverageQC.csv"

rule CODEX2_CoverageQC:
    input:
        [BAMS[sample] for sample in SAMPLES]
    output:
        Y_qc="data/work/{project}/codex2/coverageQC.csv",
        gc_qc="data/work/{project}/codex2/gc_qc.csv",
        N="data/work/{project}/codex2/library_size_factor.csv"
    params:
        project="data/work/{project}/codex2",
        bed=config['resources']['library']['targets_bed']
    script:
        "codex2_paired_getCoverage.snakemake.R"

rule CODEX2_ns:
    input:
        Y_qc="data/work/{project}/codex2/coverageQC.csv",
        gc_qc="data/work/{project}/codex2/gc_qc.csv",
        N="data/work/{project}/codex2/library_size_factor.csv"
    output:
        "data/work/{project}/codex2/chr{chr}.codex2.filtered_segments.txt"
    params:
        project="data/work/{project}/codex2",
        chr=lambda wildcards: wildcards.chr,
        normals=config['project']['normal_list']
    script:
        "codex2_paired_runChr.snakemake.R"