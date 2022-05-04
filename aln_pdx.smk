import os
import yaml
from collections import defaultdict

### FUNCTIONS ###

def map_input(wildcards):
    inputs=[]
    for RUN,_run in list(FILES[wildcards.sample].items()):
        run,lane,index=_run['PU'].split('-',2)
        inputs.append(f'bam_input/work/{wildcards.sample}/{wildcards.reference}/{run}/{lane}/{index}/mapped.bam')
    assert len(inputs)>0
    return sorted(inputs)

def get_fastqs(wildcards):
    return {'R1':'FASTQ/'+FILES[wildcards.sample][f'{wildcards.run}-{wildcards.lane}']['files'][0],'R2':'FASTQ/'+FILES[wildcards.sample][f'{wildcards.run}-{wildcards.lane}']['files'][1]}
    #fastq.yaml does not currently include file paths.

def disambiguated_fastqs(wildcards):
    return {'x':f"bam_input/work/{wildcards.sample}/bbtools/{wildcards.sample}_{config['resources']['targets_key']}_{wildcards.run}_{wildcards.lane}_{wildcards.index}.{config['reference']['key']}.fastq.gz",'y':f"bam_input/work/{wildcards.sample}/bbtools/{wildcards.sample}_{config['resources']['targets_key']}_{wildcards.run}_{wildcards.lane}_{wildcards.index}.{config['disambiguate']['key']}.fastq.gz"}

### ### PYTHON ### ###

with open(config['project']['fastq_config']) as file:
    FILES=yaml.load(file,Loader=yaml.BaseLoader)
    SAMPLES=sorted(list(FILES.keys()))
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

### ### ### RULES ### ### ###
#Consider making bams.table the desired result.

rule all:
    input:
        expand("bam_input/final/{sample}/{reference}/{sample}.bam",sample=SAMPLES,reference=config['reference']['key'])

#check if bbmap is installed?

#How do I generalize this output?
#params memory is the wrong thing. I need the actual snakemake memory
rule bbsplit:
    input:
        unpack(get_fastqs)
    output:
        x="bam_input/work/{sample}/bbtools/{sample}_{lib}_{run}_{lane}_{index}.GRCh37.fastq.gz",
        y="bam_input/work/{sample}/bbtools/{sample}_{lib}_{run}_{lane}_{index}.mm10.fastq.gz",
        refstats="metrics/{sample}/{sample}_{lib}_{run}_{lane}_{index}.refstats"
    params:
        ref_p=config['disambiguate']['ref_p'],
        memory="61440m"
    shell:
        "bbsplit.sh -Xmx{params.memory} build=1 path={params.ref_p} in1={input.R1} in2={input.R2} out_x={output.x} out_y={output.y} refstats={output.refstats}"

#input is mess, write function
rule reformat_repair:
    input:
        unpack(disambiguated_fastqs)
    output:
        R1="bam_input/work/{sample}/GRCh37/{run}/{lane}/{index}/R1.fastq.gz",
        R2="bam_input/work/{sample}/GRCh37/{run}/{lane}/{index}/R2.fastq.gz"
    params:
        memory="10240m"
    shell:
        """
        reformat.sh -Xmx{params.memory} int=t addcolon=t uniquenames=t in={input.x} out=stdout.fq |
        repair.sh -Xmx{params.memory} in=stdin.fq out1={output.R1} out2={output.R2}
        """

rule aln_pe:
    input:
        R1="bam_input/work/{sample}/GRCh37/{run}/{lane}/{index}/R1.fastq.gz",
        R2="bam_input/work/{sample}/GRCh37/{run}/{lane}/{index}/R2.fastq.gz"
    output:
        temp("bam_input/work/{sample}/GRCh37/{run}/{lane}/{index}/mapped.bam")
    params:
        LB=config['resources']['library_key'],
        fasta=config['reference']['fasta']
    threads:
        8
    shell:
        "bwa mem -M -t {threads} {params.fasta} {input.R1} {input.R2} | samtools addreplacerg -r 'ID:{wildcards.run}.{wildcards.lane}' -r 'PU:{wildcards.run}.{wildcards.lane}.{wildcards.index}' -r 'PL:illumina' -r 'LB:{params.LB}' -r 'SM:{wildcards.sample}' -@ {threads} - | samtools sort -@ {threads} -o {output}"
    #picard --metrics_accumulation_level readgroup will look at PU and will only fall back to ID if it is null'
    #Consider putting markdup in this as well. It outputs stats and multi threads. Per lane stats would be nice and HsMetrics will also do that.
    #markdup would first need samtools fixmate -m

rule input_ready:
    input:
        map_input
    output:
        temp("bam_input/work/{sample}/{reference}/sort.bam")
    params:
        input=temp("bam_input/work/{sample}/{reference}/input.bam")
    run:
        if len(input)==1:
            shell("rsync {input} {params.input}")
            shell("samtools sort {params.input} -o {output}")
        else:
            shell("samtools merge -f {params.input} {input}")
            shell("samtools sort {params.input} -o {output}")

rule MarkDuplicates:
    input:
        "bam_input/work/{sample}/{reference}/sort.bam"
    output:
        bam=temp("bam_input/work/{sample}/{reference}/mDup.bam"),
        metrics="bam_input/final/{sample}/metrics/{reference}/mark_duplicates.table"
    params:
        memory="10240m"
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/picard/2.20.7/picard.jar MarkDuplicates I={input} O={output.bam} M={output.metrics} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"

rule RealignerTargetCreator:
    input:
        "bam_input/work/{sample}/{reference}/mDup.bam"
    output:
        "bam_input/work/{sample}/{reference}/IndelRealigner.intervals"
    params:
        memory="10240m",
        reference=config['reference']['fasta'],
        known=["$HOME/resources/gatk/1000G_phase1.indels.b37.vcf","$HOME/resources/gatk/Mills_and_1000G_gold_standard.indels.b37.vcf"]
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T RealignerTargetCreator -I {input} -o {output} -known {params.known[0]} -known {params.known[1]}"

rule IndelRealigner:
    input:
        bam="bam_input/work/{sample}/{reference}/mDup.bam",
        targets="bam_input/work/{sample}/{reference}/IndelRealigner.intervals"
    output:
        temp("bam_input/work/{sample}/{reference}/realign.bam")
    params:
        memory="10240m",
        reference=config['reference']['fasta'],
        known=["$HOME/resources/gatk/1000G_phase1.indels.b37.vcf","$HOME/resources/gatk/Mills_and_1000G_gold_standard.indels.b37.vcf"]
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T IndelRealigner -I {input.bam} -o {output} -targetIntervals {input.targets} -known {params.known[0]} -known {params.known[1]}"

rule FirstPass_BaseRecalibrator:#update resources
    input:
        "bam_input/work/{sample}/{reference}/realign.bam"
    output:
        "bam_input/final/{sample}/metrics/{reference}/recal_data.table"
    params:
        memory="10240m",
        reference=config['reference']['fasta'],
        knownSites=["$HOME/resources/gatk/dbsnp_138.b37.vcf","$HOME/resources/gatk/1000G_phase1.indels.b37.vcf","$HOME/resources/gatk/Mills_and_1000G_gold_standard.indels.b37.vcf"]
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T BaseRecalibrator -I {input} -o {output} -knownSites {params.knownSites[0]} -knownSites {params.knownSites[1]} -knownSites {params.knownSites[2]}"

rule SecondPass_BaseRecalibrator:
    input:
        "bam_input/work/{sample}/{reference}/realign.bam",
        "bam_input/final/{sample}/metrics/{reference}/recal_data.table"
    output:
        "bam_input/final/{sample}/metrics/{reference}/post_recal_data.table"
    params:
        memory="10240m",
        reference=config['reference']['fasta'],
        knownSites=["$HOME/resources/gatk/dbsnp_138.b37.vcf","$HOME/resources/gatk/1000G_phase1.indels.b37.vcf","$HOME/resources/gatk/Mills_and_1000G_gold_standard.indels.b37.vcf"]
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T BaseRecalibrator -I {input[0]} -BQSR {input[1]} -o {output} -knownSites {params.knownSites[0]} -knownSites {params.knownSites[1]} -knownSites {params.knownSites[2]}"

rule AnalyzeCovariates:
    input:
        before="bam_input/final/{sample}/metrics/{reference}/recal_data.table",
        after="bam_input/final/{sample}/metrics/{reference}/post_recal_data.table"
    output:
        csv="bam_input/final/{sample}/metrics/BQSR.csv",
        pdf="bam_input/final/{sample}/metrics/BQSR.pdf"
    params:
        memory="10240m",
        reference=config['reference']['fasta'],
        knownSites=["$HOME/resources/gatk/dbsnp_135.b37.vcf","$HOME/resources/gatk/1000G_phase1.indels.b37.vcf","$HOME/resources/gatk/Mills_and_1000G_gold_standard.indels.b37.vcf"]
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T AnalyzeCovariates -before {input.before} -after {input.after} -csv {output.csv} -plots {output.pdf}"

rule PrintReads:
    input:
        bam="bam_input/work/{sample}/{reference}/realign.bam",
        bqsr="bam_input/final/{sample}/metrics/{reference}/recal_data.table"
    output:
        temp("bam_input/work/{sample}/{reference}/recal.bam")
    params:
        memory="1020m",
        reference=config['reference']['fasta']
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T PrintReads -I {input.bam} -BQSR {input.bqsr} -o {output}"

rule ValidateSamFile:
    #Errors and Warnings trigger non-zero exit
    input:
        "bam_input/work/{sample}/{reference}/recal.bam"
    output:
        "bam_input/work/{sample}/{reference}/validation_data.table"
    params:
        memory="10240m"
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/picard/2.20.7/picard.jar ValidateSamFile I={input} O={output} MODE=SUMMARY"

rule validation_pass:
    input:
        "bam_input/work/{sample}/{reference}/validation_data.table"
    output:
        "bam_input/final/{sample}/{reference}/validation_data.table"
    shell:
        """
        set +H
        if egrep -q 'No errors found' {input[0]}; then
            cp {input[0]} {output[0]}
        else
            egrep '^ERROR' {input[0]}
            exit 1
        fi
        """

rule ready_bam:
    input:
        bam="bam_input/work/{sample}/{reference}/recal.bam",
        table="bam_input/final/{sample}/{reference}/validation_data.table"
    output:
        "bam_input/final/{sample}/{reference}/{sample}.bam"
    shell:
        """
        rsync -v {input.bam} {output}
        samtools index {output}
        """
