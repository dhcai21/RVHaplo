# RVHaplo: Reconstructing Viral Haplotypes using long reads.
A viral haplotype reconstruction tool for long reads.

### Citation
Dehan Cai, Yanni Sun, Reconstructing viral haplotypes using long reads, Bioinformatics, Volume 38, Issue 8, 15 April 2022, Pages 2127–2134, https://doi.org/10.1093/bioinformatics/btac089

### E-mail: dhcai2-c@my.cityu.edu.hk
### Version: V3 (2022-08-02 updated)
##### Use [Medaka](https://github.com/nanoporetech/medaka) to polish the final result
##### Support datasets with large sizes (>100k reads) by applying MCL to subgraphs of reads. (use **-sg #_of_subgraphs**)

#### Version: V2 (2022-05-22 updated)
##### Use a C package of MCL.
##### Support multiprocessing.

## Installation
### Dependencies:
* Conda
* Python >=3.8
* samtools >= 1.4.1
* pysam
* [Medaka](https://github.com/nanoporetech/medaka)
* MCL
* Required python package: pandas >= 1.1.3 tqdm, scipy

### An easiler way to install
After cloning this respository, you can use anaconda to install the **rvhaplo.yaml** (Linux). This will install all packages you need. The command is: `conda env create -f rvhaplo.yaml -n rvhaplo`

#### An optional way to install
```
conda create -n rvhaplo python==3.8

conda activate rvhaplo

conda install -c bioconda samtools

conda install mcl

pip install medaka scipy pandas tqdm pysam
```
##### Possible problem
`'../lib/libcrypto.1.0.0.dylib' (no such file) when using samtools`

You can use the command:

`ln -s your_conda/rvhaplo/lib/libcrypto.your_exisiting_version.dylib your_conda/rvhaplo/lib/libcrypto.1.0.0.dylib`

## Usage
#### Initialization
`cd RVHaplo`<BR/>
`chmod +x rvhaplo.sh`
#### Command

`alignment: minimap2 -a reference.fasta reads.fastq > alignment.sam`

`RVHaplo:   ./rvhaplo.sh -i alignment.sam -r reference.fasta -o result -p prefix -t 8`<BR/>

`If your dataset has a large number of reads (e.g., >50k), you can use the parameter -sg n to separate the graph into n subgraphs to running.`<BR/>

`If you want to output the SNV sites only: ./rvhaplo.sh -i alignment.sam -r reference.fasta -os 1`<BR/>

`If you have multiple CPU cores, you can set the value of "-t" to accelerate the running.`<BR/>

```Usually, using the default parameters should be fine. You can check the description below to adjust the parameters. But if you don't know how to adjust the parameters, feel free to email me.`<BR/>
```

```
required arguments:
    -i  | --input:                     alignment file (sam)
    -r  | --refernece:                 reference genome (fasta)

optional arguments:
    -h  | --help :                     Print help message.
    -o  | --out :                      Path where to output the results. (default: ./result)
    -p  | --prefix STR :               Prefix of output file. (default: rvhaplo)
    -t  | --thread INT :               Number of CPU cores for multiprocessing. (default: 8)
    -mq | --map_qual INT :             Smallest mapping quality for reads. (default: 0)
    -e  | --error_rate FLOAT :         Sequencing error rate. (default: 0.1)
    -s  | --signi_level FLOAT :        Significance level for binomial tests. (default: 0.05)
    -c  | --cond_pro FLOAT :           A threshold of the maximum conditional probability for SNV sites. (default: 0.65)
    -f  | --fre_snv :                  The most dominant base' frequency at a to-be-verified site should >= fre_snv. (default: 0.80)    
    -n1 | --num_read_1 INT :           Minimum # of reads for calculating the conditional probability given one conditional site. (default: 10)
    -n2 | --num_read_2 INT :           Minimum # of reads for calculating the conditional probability given more than one conditional sites. (default: 5)
    -g  | --gap INT :                  Minimum length of gap between SNV sites for calculating the conditional probability. (default: 15)
    -ss | --smallest_snv INT :         Minimum # of SNV sites for haplotype construction. (default: 20)
    -os | --only_snv (0 or 1) :        Only output the SNV sites without running the haplotype reconstruction part. (default: 0)
    -or | --overlap_read INT :         Minimum length of overlap for creating edges between two read in the read graph. (default: 5)
    -wr | --weight_read FLOAT :        Minimum weights of edges in the read graph. (default: 0.85)
    -sg | --sub_graph INT:             Number of subgraphs to run MCL (default:1)"
    -m  | --mcl_inflation FLOAT :      Inflation of MCL algorithm. (default:2)
    -l  | --lar_cluster INT :          A threshold for seperating clusters into two groups based on sizes of clusters. (default: 50)
    -oc | --overlap_cluster INT :      A parameter related to the minimum overlap between consensus sequences of clusters. (default: 10)
    -d  | --depth INT :                Depth limitation for consensus sequences generated from clusters. (default: 5)
    -wc | --weight_cluster FLOAT :     Minimum weights between clusters in the hierarchical clustering (default: 0.8)
    -sp | --start_pos INT:             Starting position for generating consensus sequences (default: 1)
    -ep | --end_pos INT:               Ending position for generating consensus sequences. (default: 1e10)
    -a  | --abundance FLOAT :          A threshold for filtering low-abundance haplotypes. (default: 0.005)
```
`-t  | --thread`

If you have multiple CPU cores for running the tool, please use this parameter for accelerating the process. The default value is 8 (8 CPU cores).

`-mq  | --map_qual`

If you want to filter some reads with small mapping qualities, you can use this parameter (e.g., 20). When you data have a large number of reads, you can use this parameter to remove some bad-quality reads. This will help accelerate the running.

`-e  | --error_rate`

You can input a roughly estimated sequencing error rate. It will not significantly affect the SNV detection result as there are post-processing steps. And we use 0.1 as the general sequencing error rate for TGS data.

`-s  | --signi_level`

Using a small significance level value may improve the precision of detected SNV sites obtained from binomial tests. But the small significance level value may reduce recall. Thus, we suggest using the default value 0.05.

`-c  | --cond_pro`

A threshold of the maximum conditional probability for SNV sites. If the maximum conditional probability of an SNV site is less than the threshold, this site will be recognized as a fake SNV site.

`-f  | --fre_snv`

Usually, sites containing fake SNVs caused by sequencing errors still have high frequencies of the most dominant bases. And those sites with small frequencies of the most dominant base are highly possible to contain real SNVs. Thus, we only verify part of potential sites obtained from the second binomial test to accelerate the verified process.

`-n1 | --num_read_1`

Minimum number of reads for calculating the conditional probability given one conditional site. For example, P(A|B). A large value of n1 will improve the precision of detecting SNV sites, but the recall may reduce because some SNVs from low-abundant haplotypes are missed. If you are more concerned about the quality of reconstructed haplotypes instead of reconstructing the low-abundant haplotypes, you can set a large value of n1.

`-n2 | --num_read_2`

The minimum number of reads for calculating the conditional probability given more than one conditional site. For example, P(A|B1,B2,B3,...). As the number of reads covering more SNVs sites will reduce, we allow a smaller number of reads for calculating the conditional probability given more conditional sites compared to only given one conditional site. A large value of n1 will improve the precision of detecting SNV sites, but the recall may reduce because some SNVs from low-abundant haplotypes are missed. If you are more concerned about the quality of reconstructed haplotypes instead of reconstructing the low-abundant haplotypes, you can set a large value of n2. 

`-g  | --gap`

Because sites in close proximity tend to weaken the independence of sequencing errors, the distance (gap) between the target site and the given sites should be above a threshold.

`-ss | --smallest_snv`

As a small number of detected SNV sites indicates that only one strain in the sample is highly possible, the haplotype reconstruction process will stop and only output the detected SNV sites.

`-os | --only_snv`

Only output the SNV sites without running the haplotype reconstruction part.

`-or | --overlap_read`

Minimum length of overlap for creating edges between reads in the read graph. If there are many detected SNV sites and the read length is long, you can set a higher threshold. 

`-wr | --weight_read`

Minimum weights of edges in the read graph. If your data has many haplotypes or the haplotypes are highly similar, I would suggest using a higher weight, e.g., 0.9.

`-sg  | --sub_graph`

If your dataset contains more than 100k reads, the read graph will be large, which may fail the running of MCL. This paramter is to separate the read graph into subgraphs. Then, MCL will be applied to each subgraph for downsteam analysis.


`-m  | --mcl_inflation`

The parameter "Inflation" of the graph clustering algorithm Markov Cluster (MCL). Usually, using the default value 2 is enough here. For further details, please refer to https://micans.org/mcl/.

`-l  | --lar_cluster`

A threshold for seperating clusters into two groups based on sizes of clusters. Clusters with a larger size threshold will lead to more accurate consensus sequences but may miss some consensus sequences that are bridges between other consensus sequences. Users can modify the threshold based on the input number of reads. We suggest using one of (20, 30, 40, 50). 

`-oc | --overlap_cluster`

A parameter related to the minimum overlap between consensus sequences of clusters. The minimum overlap between consensus sequences is min(0.1 * x, oc), where x is the number of detected SNV sites. If there are many detected SNV sites and the read length is long, you can set a higher threshold.

`-d  | --depth`

Depth limitation for consensus sequences generated from clusters. The total number of bases (A, C, G, T) at a site should be larger than the threshold. Otherwise, ignore this site (use '-').

`-wc | --weight_cluster`

Minimum weights between clusters in the hierarchical clustering. If your data has many haplotypes or the haplotypes are highly similar, I would suggest using a higher weight, e.g., 0.9.

`-sp  | --start_pos`

Starting position for generating consensus sequences. (1-index)

`-ep  | --end_pos`

Ending position for generating consensus sequences. A large default value is for covering the whole genome. (1-index)

A threshold for filtering low-abundance haplotypes.

`-a  | --abundance`

A threshold for filtering low-abundance haplotypes.



## Output Results
All reconstructed haplotypes are summarized in a file ***_haplotypes.fasta** (Polished by Medaka). Below is an example of three haplotypes.
```
>haplotype_0_length_9181_abundance_0.50_number_of_reads_5000_depth_500
GGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGGTCTCTGGCTAACTAGGGAACC...
>haplotype_1_length_9178_abundance_0.30_number_of_reads_3000_depth_300
GTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCC...
>haplotype_2_length_9180_abundance_0.20_number_of_reads_2000_depth_200
GGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGGACC...
```

#### Clusters of Reads
For the corresponding read clusters of haplotypes, the alignment file of each cluster can be found in the folder "clusters".

#### SNV sites
You can find out the SNV sites (1-index on the reference genome) in the file "*_SNV.txt".

The base (A, C, G, T, etc.) counts from the alignments between reads and the reference genome can be found in the file "*_acgt.txt".

`
If you want to output the SNV sites only without running the haplotype reconstruction part, use the command "-os 1" or "--only_snv 1".
`

## Simulation data
Simulation datasets can be downloaded from https://drive.google.com/drive/folders/16azUqV6thGJyBThR0OaGfNgjXYhD_lEE?usp=sharing

## Tested data
The tested data can be downloaded from https://drive.google.com/drive/folders/1A4I722MO4Ujnvu_B3pugq6v6N4t1koHd?usp=sharing

`./rvhaplo.sh -i test.sam -r reference.fasta`<BR/>

