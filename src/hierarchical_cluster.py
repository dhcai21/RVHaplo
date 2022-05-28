import pickle
import numpy as np
from collections import Counter
import copy
import sys

#### functions

def start_end(c):
    res = []
    for i in c:
        temp = np.array([[Qstart[n[j]], Qend[n[j]]] for j in i])
        s = min(temp[:, 0])
        e = max(temp[:, 1])
        res.append([s, e, len(i)])
    return np.array(res)


def subcluster_seq(clus, nod, seq_mat, snv_base, major_base,depth=5):
    num_range = np.array([i for i in range(len(snv_base))])
    temp = [n[i] for i in clus]
    mat = seq_mat[temp]
    seq = np.array(['-']*len(snv_base))
    for i in range(nod[0],nod[1]+1):
        temp = Counter(mat[:, i])
        temp['-'] = 0
        max_value = max(temp.values())
        values_list = np.array(list(temp.values()))
        if sum(temp.values())>depth and sum(values_list==max_value)==1:
            nucl = max(temp, key=temp.get)
            seq[i] = nucl
    temp1 = seq!='-'
    temp2 = np.logical_or(seq == snv_base,seq==major_base)
    seq_flag = np.logical_and(temp1,temp2)
    if np.count_nonzero(seq_flag) == 0:
        seq_start_end = [num_range[-1],num_range[-1],nod[-1]]
    else:
        seq_start_end = [num_range[seq_flag][0],num_range[seq_flag][-1],nod[-1]]
    return seq, seq_flag,seq_start_end



def merge_reads(haplo, clus):
    c_all = []
    for h in haplo:
        c_temp = []
        for c in h:
            c_temp = c_temp + clus[c]
        c_all.append(c_temp)
    return c_all


def haplo_filter(haplo_fre,haplo_cluster,thres_haplo=0.005):
    fre_new = []
    clus = []
    if haplo_cluster:
        flag = haplo_fre > thres_haplo
        temp = sum(haplo_fre[flag])
        for i in range(len(flag)):
            if flag[i]:
                fre_new.append(haplo_fre[i]/temp)
                clus.append(haplo_cluster[i])
    fre_new = np.array(fre_new)
    return fre_new,clus

def weight_cluster_haplotype(seq_l, seq_m, node_l, node_m,w_thres = 0.8):
    l = len(seq_l)
    m = len(seq_m)
    index_l = np.argsort(node_l[:, 0])
    index_m = np.argsort(node_m[:, 0])
    mat = np.zeros([l, m])
    for i in range(l):
        index_i = index_l[i]
        seq_i = seq_l[index_i][0]
        seq_flag_i = seq_l[index_i][1]
        for j in range(m):
            index_j = index_m[j]
            if node_m[index_j, 0] > node_l[index_i, 1]:
                break
            seq_j = seq_m[index_j][0]
            seq_flag_j = seq_m[index_j][1]
            seq_flag_ij = np.logical_and(seq_flag_i, seq_flag_j)
            diff_flag = seq_i[seq_flag_ij] != seq_j[seq_flag_ij]
            dis = np.count_nonzero(diff_flag)
            all_l =  np.count_nonzero(seq_flag_ij)
            if all_l < overlap_num:
                continue
            temp = (1-(0.1 + dis) / (0.1+all_l))
            mat[index_i, index_j] = temp
    mat = mat * (mat>w_thres)
    return mat

def updata_haplotype(c_less, c_more, mat):
    l = len(c_less)
    c_updata = copy.deepcopy(c_more)
    for i in range(l):
        if np.count_nonzero(mat[i, :]==0) == len(c_more):
            continue
        index = np.argsort(mat[i, :])[-1]
        c_updata[index] = c_updata[index] + c_less[i]
    return c_updata


def unique_reads(clu):
    res = []
    for c in clu:
        temp = []
        for r in c:
            if r not in temp:
                temp.append(r)
        res.append(temp)
    return res

def contig_weight(seqs,nod):
    l = len(seqs)
    index = np.argsort(nod[:, 0])
    mat = np.zeros([l, l])
    for i in range(l - 1):
        index_i = index[i]
        seq_i = seqs[index_i][0]
        seq_flag_i = seqs[index_i][1]
        for j in range(i + 1, l):
            index_j = index[j]
            if nod[index_j, 0] > nod[index_i, 1]:
                break
            seq_j = seqs[index_j][0]
            seq_flag_j = seqs[index_j][1]
            seq_flag_ij = np.logical_and(seq_flag_i, seq_flag_j)
            diff_flag = seq_i[seq_flag_ij] != seq_j[seq_flag_ij]
            dis = np.count_nonzero(diff_flag)
            all_l = np.count_nonzero(seq_flag_ij)
            if all_l < overlap_num:
                continue
            temp = 1-(0.1+dis) / (0.1+all_l)
            mat[index_i,index_j] = temp
    return mat


def hierarchical_weight(clus,seq_mat,snv_base,major_base,thres = 0.8,depth = 5):
    node = start_end(clus)
    seqs = []
    node_se = []
    for i in range(len(clus)):
        temp1,temp2,temp3 = subcluster_seq(clus[i], node[i],seq_mat,snv_base,major_base,depth)
        seqs.append([temp1,temp2])
        node_se.append(temp3)
    node_se = np.array(node_se)
    mat = contig_weight(seqs,node_se)
    mat = mat*(mat>=thres)
    return mat

def hierarchical_clustering(clus,seq_mat,snv_base,major_base, thres = 0.8,depth = 5):
    num_haplotig = 1
    res_cluster = copy.deepcopy(clus)
    while True:
        contig_mat = hierarchical_weight(res_cluster,seq_mat,snv_base,major_base,thres,depth)
        maximum_weight = np.max(contig_mat)
        if maximum_weight ==0:
            break
        node_x,node_y = np.where(contig_mat == maximum_weight)
        c_merge = []
        c_merge.append(res_cluster[node_x[0]]+res_cluster[node_y[0]])
        res_cluster[node_x[0]] = []
        res_cluster[node_y[0]] = []
        res_cluster = [temp for temp in res_cluster if temp]
        for c in c_merge:
            res_cluster.append(c)
        if len(res_cluster)==num_haplotig:
            break
        else:
            num_haplotig = len(res_cluster)
    return res_cluster

####### data preparation #######
f = open(sys.argv[1], 'rb')
data = pickle.load(f)
f.close()
reads = data['reads']
Qstart = data['Qstart']
Qend = data['Qend']
seq_mat = data['seq_mat']
Coun = data['Coun']
n = np.array([i for i in range(len(Qstart))])
snv_base = data['snv_base']
major_base = data['major_base']
del data

lar_cluster = int(sys.argv[2])
depth = int(sys.argv[3])
ovlap_cluster = int(sys.argv[4])
weight_cluster = float(sys.argv[5])
fre_thres = float(sys.argv[6])
file_prefix = sys.argv[7]
f = open(f"{file_prefix}_reads_cluster.txt")
clusters = [] 
while True:
    line = f.readline().split()
    if not line:
        break
    temp = np.array([int(i) for i in line])
    clusters.append(temp)

f.close()


###  parameter initialization
_,p = seq_mat.shape
overlap_num = np.round(p*0.1)
if overlap_num > ovlap_cluster:
    overlap_num = ovlap_cluster

### seperate clusters into two groups: small & large
c_less = []
c_more = []
for c in clusters:
    if len(c) >=lar_cluster:
        c_more.append(list(c))
    else:
        c_less.append(list(c))

### hierarchical clustering
haplo_cluster = hierarchical_clustering(c_more,seq_mat,snv_base,major_base,weight_cluster,depth)

### assign small clusters to haplotype clusters
node_haplo = start_end(haplo_cluster)
node_haplo_se = []
seq_haplo = []
for i in range(len(haplo_cluster)):
    temp1,temp2,temp3 = subcluster_seq(haplo_cluster[i], node_haplo[i],seq_mat,snv_base,major_base,depth)
    seq_haplo.append([temp1,temp2])
    node_haplo_se.append(temp3)

node_haplo_se = np.array(node_haplo_se)

if len(c_less)!=0:
    node_less = start_end(c_less)
    node_se_less = []
    seq_less = []
    for i in range(len(c_less)):
        temp1, temp2, temp3 = subcluster_seq(c_less[i], node_less[i], seq_mat, snv_base, major_base, 0)
        seq_less.append([temp1, temp2])
        node_se_less.append(temp3)
    node_se_less = np.array(node_se_less)
    mat_cluster_haplotye = weight_cluster_haplotype(seq_less, seq_haplo, node_se_less, node_haplo_se, weight_cluster)
    haplo_cluster = updata_haplotype(c_less, haplo_cluster, mat_cluster_haplotye)
    node_haplo = start_end(haplo_cluster)

###  hierarchical_clustering again
haplo_cluster = hierarchical_clustering(haplo_cluster,seq_mat,snv_base,major_base,weight_cluster,depth)

###  unique reads in each cluster
haplo_cluster = unique_reads(haplo_cluster)
node_haplo = start_end(haplo_cluster)

###  abundance estimation and filter extremely low-abundant strain
haplo_fre = []
for i in node_haplo:
    haplo_fre.append(i[-1] / sum(node_haplo[:, -1]))

haplo_fre = np.array(haplo_fre)

haplo_fre,haplo_cluster = haplo_filter(haplo_fre, haplo_cluster,fre_thres)
node_haplo = start_end(haplo_cluster)


###  save result
data = {'reads': reads, 'haplo_cluster': haplo_cluster, 'haplo_fre': haplo_fre}
f = open(file_prefix+"_final.pickle", 'wb')
pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
f.close()
exit()
