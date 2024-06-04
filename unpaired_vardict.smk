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
        expand("data/work/{lib}/{sample}/vardict/variants.vep.report.csv",lib=config['resources']['targets_key'],sample=SAMPLES)

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

rule run_vep:
    input:
        "{work_dir}/{sample}/vardict/variants.vcf.gz"
    output:
        "{work_dir}/{sample}/vardict/variants.vep.vcf.gz"
    params:
        out_vcf="{work_dir}/{sample}/vardict/variants.vep.vcf",
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='/home/bwubb/.vep/Plugins/loftee',#check
        utr=config['resources']['utrannotator']
    shell:
        """
        vep -i {input} -o {params.out_vcf} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf \
        --everything \
        --canonical \
        --assembly {params.assembly} \
        --species homo_sapiens \
        --fasta {params.fa} \
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin REVEL,{params.revel} \
        --plugin SpliceAI,snv={params.splice_snv},indel={params.splice_indel} \
        --plugin gnomADc,{params.gnomAD} \
        --plugin UTRannotator,{params.utr} \
        --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \


        bgzip {params.out_vcf}
        tabix -fp vcf {output}
        """
#--plugin Downstream \
#--plugin LoF,loftee_path:{params.loftee},human_ancestor_fa:{params.human_ancestor_fa},conservation_file:{params.conservation_file},gerp_bigwig:{params.gerp_bigwig} \

rule vep_report:
    input:
        vcf="{work_dir}/{sample}/vardict/variants.vep.vcf.gz"
    output:
        csv="{work_dir}/{sample}/vardict/variants.vep.report.csv"
    shell:
        """
        python vep_vcf_parser.py -i {input.vcf} -o {output.csv} --mode single,{wildcards.sample} everything
        """
