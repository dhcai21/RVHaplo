import os
import sys
import pickle
import numpy as np
import pysam
import pandas as pd

## load clusters and initialize

file_final = sys.argv[1]
file_bam = sys.argv[2]  
file_path = sys.argv[3]   
file_acgt = sys.argv[4]
file_ref = sys.argv[5]
file_out = sys.argv[6]
s_pos = int(sys.argv[7])
e_pos = int(sys.argv[8])

attributes = ['pos', 'reads_all','deletions','A', 'C', 'G', 'T']

if file_ref == '1':
    df = pd.read_csv(file_acgt, sep='\t')
    data = df[attributes]
    mat = data.values
    chara = 'ACGT'
    flag = np.logical_and(mat[:,0]>=s_pos,mat[:,0]<=e_pos)
    mat = mat[flag,:]
    temp_depth = sum(mat[:,1])/len(mat)
    f_out = open(file_out,'w')
    bam = pysam.AlignmentFile(file_bam)
    ref_temp = bam.references[0]
    read_num = bam.count(ref_temp)
    haplotype = []
    for i in range(len(mat)):
        if mat[i,2]/mat[i,1]>=0.5:
            continue
        index = np.argsort(-mat[i,3:])[0]
        base = chara[index]
        haplotype.append(base)
    haplotype = ''.join(haplotype)
    temp_length = len(haplotype)
    f_out.write(f'>haplotype_0_length_{temp_length}_abundance_1_number_of_reads_{read_num}_depth_{temp_depth}\n')
    f_out.write(haplotype+'\n')
    f_out.close()
    exit()

f = open(file_final,'rb')
data = pickle.load(f)
f.close()
reads = data['reads']
n = np.array([i for i in range(len(reads))])
haplo_cluster = data['haplo_cluster']
haplo_fre = np.array(data['haplo_fre'])
reads = np.array(reads)
index = list(np.argsort(-haplo_fre))
haplo_fre = haplo_fre[index]
haplo_cluster = [haplo_cluster[i] for i in index]
del data
def out_reads(reads,n,haplo):
    R = []
    for h in haplo:
        r_temp = reads[list(n[h])]
        R.append(r_temp)
    return R

R = out_reads(reads,n,haplo_cluster)


if len(haplo_fre)==1:
    f_out = open(file_out,'w')
    df = pd.read_csv(file_acgt, sep='\t')
    data = df[attributes]
    mat = data.values
    chara = 'ACGT'
    flag = np.logical_and(mat[:,0]>=s_pos,mat[:,0]<=e_pos)
    mat = mat[flag,:]
    temp_depth = sum(mat[:,1])/len(mat)
    bam = pysam.AlignmentFile(file_bam)
    ref_temp = bam.references[0]
    read_num = bam.count(ref_temp)
    haplotype = []
    for i in range(len(mat)):
        if mat[i,2]/mat[i,1]>=0.5:
            continue
        index = np.argsort(-mat[i,3:])[0]
        base = chara[index]
        haplotype.append(base)
    haplotype = ''.join(haplotype)
    temp_length = len(haplotype)
    f_out.write(f'>haplotype_0_length_{temp_length}_abundance_1_number_of_reads_{read_num}_depth_{temp_depth}\n')
    f_out.write(haplotype+'\n')
    f_out.close()
    exit()

f_out = open(file_out,'w')
for i in range(len(R)):
    bam = pysam.AlignmentFile(file_bam)
    qnames = R[i]
    read_num = 0
    out_name = file_path+"/clusters/cluster_"+str(i)+".bam"
    obam = pysam.AlignmentFile(out_name, "wb", template= bam)
    for b in bam.fetch(until_eof=True):
        if b.query_name in qnames:
            obam.write(b)
            read_num = read_num + 1
    obam.close()
    bam.close()
    file_cluster_bam = out_name
    file_cluster_bam_sorted = file_cluster_bam[0:-4]+'_sorted.bam'
    os.system(f"samtools sort {file_cluster_bam} -o {file_cluster_bam_sorted}")
    os.system(f"samtools index {file_cluster_bam_sorted}")
    os.system(f"python ./src/count_frequency.py {file_cluster_bam_sorted} {file_path}/clusters/cluster_{i}_acgt.txt")
    seq_temp = []
    file_acgt = file_path+"/clusters/cluster_"+str(i) + '_acgt.txt'
    df = pd.read_csv(file_acgt, sep='\t')
    data = df[attributes]
    mat = data.values
    flag = np.logical_and(mat[:,0]>=s_pos,mat[:,0]<=e_pos)
    mat = mat[flag,:]
    temp_depth = sum(mat[:,1])/len(mat)
    for j in range(len(mat)):
        pos = mat[j, 0] - 1
        if mat[j,2]/mat[j,1]>=0.5:
            continue
        index = np.argsort(-mat[j,3:])[0]
        base = chara[index]
        seq_temp.append(base)
    h = ''.join(seq_temp)
    temp_length = len(h)
    f_out.write(f'>haplotype_{i}_length_{temp_length}_abundance_{haplo_fre[i]}_number_of_reads_{read_num}_depth_{temp_depth}\n')
    f_out.write(h+'\n')

f_out.close()
exit()
