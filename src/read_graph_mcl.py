import sys
import pysam
import numpy as np
import pickle
import networkx as nx
import markov_clustering as mc

### load data
file_in=sys.argv[1]
file_snv = sys.argv[2]
cond_pro = float(sys.argv[3])
smallest_snv = int(sys.argv[4])
num_read_1 = int(sys.argv[5])
num_read_2 = int(sys.argv[6])
gap = int(sys.argv[7])
weight_read = float(sys.argv[8])
ovlap_read = int(sys.argv[9])
mcl_inflation = float(sys.argv[10])
file_prefix = sys.argv[11]


bamfile = pysam.AlignmentFile(file_in)
ref_name = bamfile.references[0]
f_snv = open(file_snv,'r')
snv_sites = f_snv.readline().split()
snv_sites = np.array([int(i)-1 for i in snv_sites])
f_snv.close()


### data initialization
def error_rate(qual):
    er = 10 ** (-qual / 10)
    return er

Base = {}
Coun = []
Map_acc = {}
Reads_all = []
Qstart = []
Qend = {}
s = 0
for column in bamfile.pileup(ref_name, min(snv_sites), max(snv_sites)+1, truncate=True, min_base_quality=0):
    if column.pos in snv_sites:
        i = 0
        mq = column.get_mapping_qualities()
        coun = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for reads in column.pileups:
            read_name = reads.alignment.query_name
            pos = reads.query_position
            Map_acc[read_name] = 1-error_rate(mq[i])
            if reads.is_del or reads.is_refskip:
                base_temp = '-'
            else:
                base_temp = reads.alignment.query_sequence[pos]
            if read_name not in Reads_all:
                Reads_all.append(read_name)
                Base[read_name] = {column.pos:base_temp}
                Qstart.append(s)
                Qend[read_name] = s
            else:
                Base[read_name][column.pos] = base_temp
                Qend[read_name] = s
            if base_temp in list(coun.keys()):
                coun[base_temp] = coun[base_temp] + 1
            i = i + 1
        Coun.append(list(coun.values()))
        s = s + 1

Coun = np.array(Coun)
Qstart = np.array(Qstart)

###  generate the sequence matrix
seq_mat = []
num_snv = len(snv_sites)
num_reads = len(Reads_all)
for i in range(num_reads):
    read_name = Reads_all[i]
    temp = ['-']*Qstart[i]+list(Base[read_name].values())+['-']*(num_snv - Qend[read_name]-1)
    seq_mat.append(temp)

seq_mat = np.array(seq_mat)
del Base

####  calculate nucleotide frequency and rank
def rank_nucl(Coun):
	res = []
	for i in Coun:
		res.append(list(np.sort(i)/sum(i)))
	return res

nucl_fre = np.array(rank_nucl(Coun))

###   find out the candicate site to be verified >=0.8  
to_be_verified_sites = np.array([i for i in range(num_snv)])
to_be_verified_sites = to_be_verified_sites[nucl_fre[:,-1]>=0.8]

###  calculate the 1-order  probability for the SNV at each candicate site
subsitution_rate = {i:nucl_fre[i,2]/sum(nucl_fre[i,2:]) for i in to_be_verified_sites}

###  reserve the top-2 bases  
chara = ['A','C','G','T']
snv_base = []
major_base = []
for i in range(len(Coun)):
    temp = np.argsort(Coun[i])
    snv_base.append(chara[temp[-2]])
    major_base.append(chara[temp[-1]])

snv_base = np.array(snv_base)
major_base = np.array(major_base)

seq_snv_flag = []
seq_major_flag = []
for i in range(len(seq_mat)):
    temp = seq_mat[i]==snv_base
    seq_snv_flag.append(temp)
    temp = seq_mat[i]==major_base
    seq_major_flag.append(temp)

seq_snv_flag = np.array(seq_snv_flag)
seq_major_flag = np.array(seq_major_flag)

###   calculate the conditional probability
def cond_pro_1_order(c_site,g_site,ssf,smf,v_min = 10):
    #c_site: current site
    #g_site: given site
    #ssf: seq_snv_flag
    #smf: seq_major_flag
    current_snv_flag = ssf[:,c_site]
    current_major_flag = smf[:,c_site]
    given_flag = ssf[:,g_site]
    num_snv = np.count_nonzero(current_snv_flag[given_flag])
    num_major = np.count_nonzero(current_major_flag[given_flag])
    denominator  = num_snv + num_major
    numerator = num_snv
    if denominator > v_min:
        pro = numerator/(denominator+1e-16)
    else:
        pro = 0   
    return pro

def cond_pro_extend(c_site,g_site,cg_flag,ssf,smf):
    #c_site: current site
    #g_site: given site
    #ssf: seq_snv_flag
    #smf: seq_major_flag
    #current_given_flag
    current_snv_flag = ssf[:,c_site]
    current_major_flag = smf[:,c_site]
    given_flag = np.logical_and(ssf[:,g_site],cg_flag)
    num_snv = np.count_nonzero(current_snv_flag[given_flag])
    num_major = np.count_nonzero(current_major_flag[given_flag])
    denominator = num_snv + num_major
    numerator = num_snv
    pro = numerator/(denominator+1e-16)
    return pro,given_flag,denominator

def verify_snv_site(to_be_verified_sites,seq_snv_flag,seq_major_flag,v_min_1=10,v_min_2=5,hd=15):
    snv_pro = []
    for i in to_be_verified_sites:
        temp_flag = to_be_verified_sites!=i
        con_pro = []
        given_sites = to_be_verified_sites[temp_flag]
        for j in given_sites:
            pro = cond_pro_1_order(i,j,seq_snv_flag,seq_major_flag,v_min_1)
            con_pro.append(pro)
        con_pro = np.array(con_pro)
        index_sort = np.argsort(-con_pro)
        current_pro = subsitution_rate[i]
        current_given_flag = np.array([True]*num_reads)
        for k in index_sort:
            given_site = given_sites[k]
            if abs(snv_sites[given_site]-snv_sites[i]) < hd:
                continue
            new_pro,new_given_flag,new_v = cond_pro_extend(i,given_site,current_given_flag,seq_snv_flag,seq_major_flag)
            if new_pro>current_pro and new_v>=v_min_2:
                current_pro = new_pro
                current_given_flag = new_given_flag
            else:
                break
        snv_pro.append(current_pro)
    return snv_pro

snv_max_pro = verify_snv_site(to_be_verified_sites,seq_snv_flag,seq_major_flag,num_read_1,num_read_2,gap)
snv_max_pro = np.array(snv_max_pro)
fake_snv_sites = to_be_verified_sites[snv_max_pro<cond_pro]
final_snv = [i for i in range(num_snv) if i not in fake_snv_sites]

###  update snv_sites file
f = open(file_snv,'w')
for i in final_snv:
    f.write(str(snv_sites[i]+1)+'\t')


### check the satisfaction of SNV sites's number
if len(final_snv)<smallest_snv:
    print("There is no detected SNV site.\nexit")
    f.write("\nZero SNV sites indicates one haplotype\nexit")
    f.close()
    exit()

###  update some variables
Coun = Coun[final_snv]
seq_mat = seq_mat[:,final_snv]
snv_base = snv_base[final_snv]
major_base = major_base[final_snv]
Qstart = []
Qend = []
seq_flag = []
snv_range = np.array([i for i in range(len(final_snv))])
for i in range(len(seq_mat)):
    temp1 = seq_mat[i]!='-'
    temp2 = np.logical_or(seq_mat[i]==snv_base,seq_mat[i]==major_base)
    if np.count_nonzero(temp1)==0:
        Qstart.append(snv_range[-1])
        Qend.append(snv_range[-1])
    else:
        Qstart.append(snv_range[temp1][0])
        Qend.append(snv_range[temp1][-1])
    seq_flag.append(np.logical_and(temp1,temp2))

Qstart = np.array(Qstart)
Qend = np.array(Qend)
seq_flag = np.array(seq_flag)
index_r = np.argsort(Qstart)

print("starting constructing the read graph")

def cd_score(Reads_all, Map_acc,Qstart,Qend,seq_mat,seq_flag,index_r,snv_overlap=5, threshold=0.8):
    W = nx.Graph()
    m = len(index_r)
    reads = []
    flag_ij = []
    for i in range(m-1):
        index_i = index_r[i]
        read_name_i = Reads_all[index_i]
        reads.append(read_name_i)
        qstart_i = Qstart[index_i]
        qend_i = Qend[index_i]
        qmid_i = int((qstart_i + qend_i)/2)
        seq_i = seq_mat[index_i]
        seq_flag_i = seq_flag[index_i]
        if int(i/1000)>=0 and (i%1000)==0:
            print(i,m)
        for j in range(i+1,m):
            index_j = index_r[j]
            qend_j = Qend[index_j]
            qstart_j = Qstart[index_j]
            read_name_j = Reads_all[index_j]
            if qstart_j > qmid_i:
                break
            seq_j = seq_mat[index_j]
            seq_flag_j = seq_flag[index_j]
            seq_flag_ij = np.logical_and(seq_flag_i, seq_flag_j)
            diff_flag = (seq_i[seq_flag_ij]!=seq_j[seq_flag_ij])
            d_score = np.count_nonzero(diff_flag)
            N = np.count_nonzero(seq_flag_ij)
            p_cor_i = Map_acc[read_name_i]
            p_cor_j = Map_acc[read_name_j]
            if N >= snv_overlap:
                weight_temp = (p_cor_i*p_cor_j)*(1-(0.1+d_score)/(N+0.1))
                if weight_temp >= threshold:
                    W.add_edge(i,j,weight=weight_temp)
    return W,reads

G,reads = cd_score(Reads_all, Map_acc,Qstart,Qend,seq_mat,seq_flag,index_r,ovlap_read,weight_read)
reads.append(Reads_all[index_r[-1]])
Qstart = Qstart[index_r]
Qend = Qend[index_r]
seq_mat = seq_mat[index_r]
n = np.array(G.nodes())
matrix = nx.to_scipy_sparse_matrix(G)

### MCL clustering
result = mc.run_mcl(matrix,inflation=mcl_inflation)
clusters = mc.get_clusters(result)
del matrix
data = {'seq_mat':seq_mat,'n':n,'reads':reads,'Qstart':Qstart,'Qend':Qend,'Clusters':clusters,'Coun':Coun,'snv_base':snv_base,'major_base':major_base}
f_pic = open(file_prefix+'_graph.pickle','wb')
pickle.dump(data,f_pic,protocol=pickle.HIGHEST_PROTOCOL)
f_pic.close()
del data
del G
exit()
