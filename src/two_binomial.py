import sys
import pandas as pd
from tqdm import tqdm
import concurrent.futures
import numpy as np
from scipy import stats
import numpy as np

#### input parameters
error = float(sys.argv[1])
threshold = float(sys.argv[2])
file_acgt = sys.argv[3]
f_out = open(sys.argv[4],'w')
thread = int(sys.argv[5])
s_pos = int(sys.argv[6])
e_pos = int(sys.argv[7])

### load data
df = pd.read_csv(file_acgt,sep='\t')
attributes = ['pos','reads_all','deletions','A','C','G','T']
data = df[attributes]
mat = data.values
flag = np.logical_and(mat[:,0]>=s_pos,mat[:,0]<=e_pos)
mat = mat[flag,:]
pos_length = len(mat)

### first binomal test
def first_binomial(i,mat=mat,error=error):
    x_value = max(mat[i,2:])
    n_value = mat[i,1]
    p_value = stats.binom_test(x_value,n=n_value,p=1-error,alternative='less')
    return p_value

query_list = [i for i in range(pos_length)]
print("first binomial test")
with concurrent.futures.ProcessPoolExecutor(thread) as executor:
    p_values = list(tqdm(executor.map(first_binomial, query_list), total=len(query_list)))


normal_sites = []
for i in range(pos_length):
    if p_values[i] >= threshold:
        normal_sites.append(i)

### estimate beta
total = []
weight = []
for site in normal_sites:
    if sum(mat[site,3:]) == 0:
        continue
    temp = mat[site,3:]/sum(mat[site,3:])
    temp = temp[np.argsort(-temp)]
    weight.append(temp)
    total.append(mat[site,1])
weight = np.array(weight)
tol = np.array(total)/sum(total)
for i in range(len(tol)):
    weight[i,:] = weight[i,:]*tol[i]
ref_init = np.sum(weight,0)
beta = ref_init[1]

### second binomial test
def second_binomial(i,mat=mat,beta=beta):
    temp = mat[i,3:]
    obes = temp[np.argsort(-temp)]
    x_value = obes[1]
    n_value = np.sum(obes)
    p_value = stats.binom_test(x_value,n=n_value,p=beta,alternative='greater')
    return p_value

print("second binomial test")
with concurrent.futures.ProcessPoolExecutor(thread) as executor:
    p_values = list(tqdm(executor.map(second_binomial, query_list), total=len(query_list)))

SNV_final = []
for i in range(pos_length):
    if p_values[i] < threshold:
        SNV_final.append(mat[i,0])
        f_out.write(str(mat[i,0])+'\t')

if not SNV_final:
    print("There is no detected SNV site.\nexit")
    f_out.write("\nZero SNV sites indicates one haplotype\nexit")

f_out.close()
exit()
