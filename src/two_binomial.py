import sys
import pandas as pd
import numpy as np
from scipy import stats
import numpy as np

#### input parameters
error = float(sys.argv[1])
threshold = float(sys.argv[2])
file_acgt = sys.argv[3]
f_out = open(sys.argv[4],'w')

### load data
df = pd.read_csv(file_acgt,sep='\t')
attributes = ['pos','reads_all','deletions','A','C','G','T']
data = df[attributes]
mat = data.values
pos_length,_ = data.shape

### first binomal test
normal_sites = []
for i in range(pos_length):
    x_value = max(mat[i,2:])
    n_value = mat[i,1]
    p_value = stats.binom_test(x_value,n=n_value,p=1-error,alternative='less')
    if p_value >= threshold:
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
beta = ref_init[i]

### second binomial test
SNV_final = []
for i in range(pos_length):
    temp = mat[i,3:]
    obes = temp[np.argsort(-temp)]
    x_value = obes[1]
    n_value = np.sum(obes)
    p_value = stats.binom_test(x_value,n=n_value,p=beta,alternative='greater')
    if p_value < threshold:
        SNV_final.append(mat[i,0])
        f_out.write(str(mat[i,0])+'\t')

if not SNV_final:
    print("There is no detected SNV site.\nexit")
    f_out.write("\nZero SNV sites indicates one haplotype\nexit")

f_out.close()
exit()
