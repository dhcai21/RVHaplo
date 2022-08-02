import os
import sys
import numpy as np
from collections import Counter
import pickle

file_path = sys.argv[1]
file_prefix = sys.argv[2]

os.system(f'rm -rf {file_path}/medaka/fastx')
os.system(f'mkdir -p {file_path}/medaka/fastx')
f = open(f'{file_path}/{file_prefix}_clusters.pickle','rb')
temp = Counter(pickle.load(f)['haplo_fre'])
num = len(temp)
f.close()
f = open(f'{file_path}/{file_prefix}_consensus.fasta','r')
for i in range(num):
    os.system(f'samtools fastq {file_path}/clusters/cluster_{i}.bam > {file_path}/medaka/fastx/cluster_{i}.fastq')
    f_o = open(f'{file_path}/medaka/fastx/consensus_{i}.fasta','w')
    f_o.write(f.readline())
    f_o.write(f.readline())
    f_o.close()

exit()
        
    
