import yaml
import re
import glob
import os
from collections import defaultdict


def infer_normal_cohorts():
    normals_cohort_input=defaultdict(set)
    fastq=re.compile('([a-zA-Z0-9-]+)_S0736811_(FGC\d{4})_([1-8])_[ACGT-]+_R1\.fastq\.gz')
    all_normals=glob.glob('FASTQ/KT*.fastq.gz')
    for file in all_normals:
        match=fastq.match(os.path.basename(file))
        if match:
            #print(match.group(1,2,3))
            normals_cohort_input['.'.join(list(match.group(2,3)))].add(f"data/work/{config['resources']['targets_key']}/{match.group(1)}/mutect2/normal.vcf.gz")
    return normals_cohort_input


def infer_cohort_match():
    #matched_cohort=defaultdict(str)
    matched_cohort={}
    fastq=re.compile('([a-zA-Z0-9-]+)_S0736811_(FGC\d{4})_([1-8])_[ACGT-]+_R1\.fastq\.gz')#HERE
    if os.path.exists('tumors.list'):
        not_normals=[]
        with open('tumors.list','r') as file:
            lines=file.read().splitlines()
        for i in lines:
            not_normals+=glob.glob(f'FASTQ/{i}_*.fastq.gz')
    else:
        not_normals=glob.glob('FASTQ/WM*.fastq.gz')
        #if empty raise
    for file in not_normals:
        match=fastq.match(os.path.basename(file))
        if match:
            #print(match.group(1,2,3))
            assert match.group(1) not in matched_cohort
            matched_cohort[match.group(1)]=f"data/work/{config['resources']['targets_key']}/{match.group(2)}.{match.group(3)}/mutect2/pon.vcf.gz"
    #for k,v in matched_cohort.items():
    #    print(k,v)
    return matched_cohort

#if config['project']['pon_config']:
#    pass
#    #do something
#else:
#    with open(config['project']['fastq_config'],'r') as file:
#        FASTQ=yaml.load(file)


#normals_cohort_input=infer_normal_cohorts()
#assert len(normals_cohort_input)>0
#matched_cohort=infer_cohort_match()
with open('tumor.list','r') as file:
    TUMORS=file.read().splitlines()

def bam_input(wildcards):
    try:
        return 'bam_input/final/{sample}/{reference}/{sample}.ready.bam'.format(sample=wildcards.sample,reference=config['reference']['key'])
    except AttributeError:
        return 'bam_input/final/{tumor}/{reference}/{tumor}.ready.bam'.format(tumor=wildcards.tumor,reference=config['reference']['key'])


#def pon_input(wildcards):
#    return normals_cohort_input[wildcards.cohort]

def pon_input(wildcards):
    #with open(config['project']['pon_config'],'r') as p:
    #    PON=yaml.load(p,Loader=yaml.BaseLoader)
    with open(f"{config['project']['normal_list']}",'r') as file:
        NORMALS=file.read().splitlines()
    return [f"data/work/{config['resources']['targets_key']}/{sample}/mutect2/normal.vcf.gz" for sample in NORMALS]


with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def sample_bam(wildcards):
    return BAMS[wildcards.sample]

def tumor_bam(wildcards):
    return BAMS[wildcards.tumor]

def unpaired_input(wildcards):
    return matched_cohort[wildcards.tumor]

wildcard_constraints:
    work_dir=f"data/work/{config['resources']['targets_key']}",
    results_dir=f"data/final/{config['project']['name']}"
#I was unable to run initial_mutect2_results() when +COSMIC was included in the project name

rule all:
    input:
        expand('data/work/{targets}/{tumor}/annovar/mutect2.hg19_multianno.report.tsv',tumor=TUMORS,targets=config['resources']['targets_key'])

rule normal_detection:
    input:
        sample_bam
    output:
        '{work_dir}/{sample}/mutect2/normal.vcf.gz'
    params:
        reference=config['reference']['fasta'],
        intervals=config['resources']['targets_intervals']
    #cosmic file
    shell:
        'java -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T MuTect2 -I:tumor {input} --artifact_detection_mode -L {params.intervals} -o {output}'

#make list that shares with codex2?
rule pon_creation:
    input:
        pon_input
    output:
        'project_data/mutect2/pon.vcf.gz'
    params:
        V=lambda wildcards,input:' '.join(['--variant %s' % v for v in input]),
        reference=config['reference']['fasta'],
        intervals=config['resources']['targets_intervals']
    shell:
        'java -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T CombineVariants {params.V} -minN 2 --setKey "null" --filteredAreUncalled --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED -L {params.intervals} -o {output}'

rule Unpaired_Mutect2:
    input:
        tumor=tumor_bam,
        pon='project_data/mutect2/pon.vcf.gz'
    output:
        '{work_dir}/{tumor}/mutect2/somatic.vcf.gz'
    params:
        reference=config['reference']['fasta'],
        intervals=config['resources']['targets_intervals']
    shell:
        'java -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T MuTect2 -I:tumor {input.tumor} --normal_panel {input.pon} -L {params.intervals} -o {output}'

    #add --activeRegionOut Output the active region to this IGV formatted file
    #add --bamOutput

rule Select_Filter:
    input:
        '{work_dir}/{tumor}/mutect2/somatic.vcf.gz'
    output:
        '{work_dir}/{tumor}/mutect2/somatic.pass.vcf.gz'
    params:
        reference=config['reference']['fasta'],
        intervals=config['resources']['targets_intervals']
    shell:
        #aU?
        'bcftools view -f PASS -O z -o {output} {input}'

rule annotate_mutect2:
    input:
        "{work_dir}/{tumor}/mutect2/somatic.pass.vcf.gz"
    output:
        "{work_dir}/{tumor}/annovar/mutect2.hg19_multianno.vcf.gz"
    params:
        out_p="{work_dir}/{tumor}/annovar/mutect2",
        humandb="/home/bwubb/resources/annovar/humandb"
    shell:
        """
        table_annovar.pl {input} {params.humandb} --buildver hg19 --vcfinput --outfile {params.out_p} --protocol refGene,cytoband,genomicSuperDups,dbscsnv11,avsnp150,dbnsfp35a,mcap,revel,popfreq_max_20150413,exac03,exac03nontcga,gnomad211_exome,intervar_20180118,icgc21,cosmic84_coding,cosmic84_noncoding,clinvar_20190305 --operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -remove
        bgzip {params.out_p}.hg19_multianno.vcf
        tabix -fp vcf {output}
        """

rule run_annovartools_mutect2:
    input:
        "{work_dir}/{tumor}/annovar/mutect2.hg19_multianno.vcf.gz"
    output:
        "{work_dir}/{tumor}/annovar/mutect2.hg19_multianno.report.tsv"
    params:
        header='/home/bwubb/resources/annovar/annotation-header.20210208.txt',
        mode='single'
    shell:
        'python annovartools.v04.py -I {input} -O {output} --header {params.header} -m {params.mode} {wildcards.tumor}'