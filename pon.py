import glob
import re
import os
from collections import defaultdict

cohorts=defaultdict(set)
fastq=re.compile('([a-zA-Z0-9-]+)_(FGC\d{4})_([1-8])_1_[ACGT-]+\.fastq\.gz')
all_normals=glob.glob('FASTQ/KT*.fastq.gz')

#def determine_cohorts:
for file in all_normals:
    match=fastq.match(os.path.basename(file))
    if match:
        #print(match.group(1,2,3))
        cohorts['.'.join(list(match.group(2,3)))].add(match.group(1))

#def cohort_list:
for k,v in cohorts.items():
    print(k,sorted(list(v)))
    with open('{}.list'.format(k),'w') as file:
        file.write('\n'.join(['data/final/{}/MelanomaTargeted_v3.S03094091/gatk-haplotype.g.vcf.gz'.format(x) for x in sorted(list(v))]))