import sys
import pysam
import numpy as np

### load data
file_in=sys.argv[1]
file_out = sys.argv[2]

bamfile = pysam.AlignmentFile(file_in)
ref_name = bamfile.references[0]
length = bamfile.lengths[0]
a,c,g,t = bamfile.count_coverage(ref_name,quality_threshold = 0)


bamfile = pysam.AlignmentFile(file_in)
reads_all = []
pos_all = []
for column in bamfile.pileup(ref_name,0,length,min_base_quality=0,max_depth=5000000):
    reads_all.append(column.n)
    pos_all.append(column.pos)

reads_all = np.array(reads_all)
pos_all = np.array(pos_all)


f = open(file_out,'w')
f.write("chrom\tpos\treads_all\tdeletions\tA\tC\tG\tT\n")
for i in range(len(pos_all)):
    index = pos_all[i]
    depth = reads_all[i]
    pos = index + 1
    if depth>0:
        deletion = depth-a[index]-c[index]-g[index]-t[index]
        f.write(f"{ref_name}\t{pos}\t{depth}\t{deletion}\t{a[index]}\t{c[index]}\t{g[index]}\t{t[index]}\n")

f.close()
exit()

