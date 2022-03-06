import os

try:
    s=open(config.get('project',{}).get('tumor_list','tumor.list'),'r')
except FileNotFoundError:
    s=open(config.get('project',{}).get('sample_list','sample.list'),'r')
SAMPLES=s.read().splitlines()
s.close()
for sample in SAMPLES:
    os.makedirs(f'logs/cluster/{sample}',exist_ok=True)


with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def sample_bam(wildcards):
    return BAMS[wildcards.sample]

rule all:
    input:
        expand("data/work/{lib}/{sample}/annovar/vardict.hg19_multianno.report.tsv",lib=config['resources']['targets_key'],sample=SAMPLES)

rule run_VarDictJava:#First a more lenient -P val, not sure what
    input:
        sample_bam
    output:
        "{work_dir}/{sample}/vardict/variants.vcf.gz"
    params:
        init="{work_dir}/{sample}/vardict/raw.vcf.gz",
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed'],
        path="/home/bwubb/software/VarDictJava/VarDict",
        AF_THR=0.01
    shell:
        """
        /home/bwubb/software/VarDictJava/build/install/VarDict/bin/VarDict -G {params.ref} -f {params.AF_THR} -N {wildcards.sample} -b {input} -c 1 -S 2 -E 3 -g 4 {params.bed} | {params.path}/teststrandbias.R | {params.path}/var2vcf_valid.pl -N {wildcards.sample} -f {params.AF_THR} | bgzip -c > {params.init}
        tabix -fp vcf {params.init}
        gatk UpdateVCFSequenceDictionary -V {params.init} --source-dictionary {input} --output {output}
        """

rule annotate_vardict:
    input:
        "{work_dir}/{sample}/vardict/variants.vcf.gz"
    output:
        "{work_dir}/{sample}/annovar/vardict.hg19_multianno.vcf.gz"
    params:
        out_p="{work_dir}/{sample}/annovar/vardict",
        humandb="/home/bwubb/resources/annovar/humandb"
    shell:
        """
        table_annovar.pl {input} {params.humandb} --buildver hg19 --vcfinput --outfile {params.out_p} --protocol refGene,cytoband,genomicSuperDups,dbscsnv11,avsnp150,dbnsfp35a,mcap,revel,popfreq_max_20150413,exac03,exac03nontcga,gnomad211_exome,intervar_20180118,icgc21,cosmic84_coding,cosmic84_noncoding,clinvar_20190305 --operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -remove
        bgzip {params.out_p}.hg19_multianno.vcf
        tabix -fp vcf {output}
        """

rule run_annovartools_vardict:
    input:
        "{work_dir}/{sample}/annovar/vardict.hg19_multianno.vcf.gz"
    output:
        "{work_dir}/{sample}/annovar/vardict.hg19_multianno.report.tsv"
    params:
        header='/home/bwubb/resources/annovar/annotation-header.20210208.txt',
        mode='single'
    shell:
        'python annovartools.v03.py -I {input} -O {output} --header {params.header} -m {params.mode} {wildcards.sample}'