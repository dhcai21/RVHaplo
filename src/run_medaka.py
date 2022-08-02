import os
import sys
import numpy as np
from collections import Counter
import pickle

file_path = sys.argv[1]
file_prefix = sys.argv[2]


f = open(f'{file_path}/{file_prefix}_clusters.pickle','rb')
temp = Counter(pickle.load(f)['haplo_fre'])
num = len(temp)
f.close()
for i in range(num):
    os.system(f"rm -rf {file_path}/medaka/medaka_{i}")
    os.system(f'medaka_consensus -i {file_path}/medaka/fastx/cluster_{i}.fastq -d {file_path}/medaka/fastx/consensus_{i}.fasta -o {file_path}/medaka/medaka_{i}')

os.system(f"cat {file_path}/medaka/medaka_{0}/consensus.fasta > {file_path}/{file_prefix}_haplotypes.fasta")
for i in range(1,num):
    os.system(f"cat {file_path}/medaka/medaka_{i}/consensus.fasta >> {file_path}/{file_prefix}_haplotypes.fasta")

exit()
