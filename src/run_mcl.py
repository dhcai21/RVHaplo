import os
import sys
file_prefix = sys.argv[1]
thread = sys.argv[2]
mcl_inflation = sys.argv[3]
sub_graph = sys.argv[4]

def mcl(file_prefix=file_prefix,thread=thread,mcl_inflaction=mcl_inflation):
    os.system(f'mcxload -abc {file_prefix}_reads_graph.txt --stream-mirror --write-binary -o {file_prefix}_reads_graph.mci -write-tab {file_prefix}_reads_graph.tab')
    os.system(f'rm {file_prefix}_reads_graph.txt')
    os.system(f'mcl {file_prefix}_reads_graph.mci -te {thread} -I {mcl_inflation} -l 1 -L 100 -o {file_prefix}_mcl_result.icl')
    os.system(f'rm {file_prefix}_reads_graph.mci')
    os.system(f'mcxdump -icl {file_prefix}_mcl_result.icl -o {file_prefix}_reads_cluster.txt -tabr {file_prefix}_reads_graph.tab')
    os.system(f'rm {file_prefix}_mcl_result.icl')
    os.system(f'rm {file_prefix}_reads_graph.tab')

if sub_graph == '1':
    mcl(file_prefix=file_prefix)
else:
    for i in range(int(sub_graph)):
        mcl(file_prefix=f'{file_prefix}{i}')
    os.system(f'cat {file_prefix}0_reads_cluster.txt > {file_prefix}_reads_cluster.txt')
    os.system(f'rm {file_prefix}0_reads_cluster.txt')
    for i in range(1,int(sub_graph)):
        os.system(f'cat {file_prefix}{i}_reads_cluster.txt >> {file_prefix}_reads_cluster.txt')
        os.system(f'rm {file_prefix}{i}_reads_cluster.txt')



