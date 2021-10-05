#! /bin/bash
### required arguments
file_sam=""
file_ref=""
### optional arguments
file_path='./'
prefix="rvhaplo"
error_rate=0.1
signi_level=0.05
cond_pro=0.65
fre_snv=0.8
num_read_1=10
num_read_2=5
gap=15
smallest_snv=20
ovlap_read=5
weight_read=0.8
mcl_inflation=2
lar_cluster=50
ovlap_cluster=10
depth=5
weight_cluster=0.8
abundance=0.005


function help_info() {
	echo "Usage: $0 -i alignment.sam -r ref_genome.fasta [options]"
	echo ""
	echo "RVHaplo: Reconstructing viral haplotypes using long reads"
	echo ""
	echo "Author: Dehan CAI"
	echo "Date:   Sept 2021"
	echo ""
	echo "    -i | --input:                     alignment file (sam)"
	echo "    -r | --refernece:                 reference genome (fasta)"
	echo ""
	echo "    Options:"
	echo "    -o  | --out:                      Path where to output the results. (default:./)"
	echo "    -p  | --prefix STR:               Prefix of output file. (default: rvhaplo)"
	echo "    -e  | --error_rate FLOAT:         Sequencing error rate. (default: 0.1)"
	echo "    -s  | --signi_level FLOAT:        Significance level for binomial tests. (default: 0.05)"
	echo "    -c  | --cond_pro FLOAT:           A threshold for the maximum conditional probability of a SNV site. (default: 0.65)"
	echo "    -f  | --fre_snv FLOAT:            The most dominant base' frequency at a to-be-verified site should >= fre_snv. (default: 0.80)"
	echo "    -n1 | --num_read_1 INT:           Minimum # of reads for calculating the conditional probability given one conditional site. (default:10)"
	echo "    -n2 | --num_read_2 INT:           Minimum # of reads for calculating the conditional probability given more than one conditional sites. (default: 5)"
	echo "    -g  | --gap INT:                  Minimum length of gap between SNV sites for calculating the conditional probability. (default:15)"
	echo "    -s  | --smallest_snv INT:         Minimum # of SNV sites for haplotype construction. (default:20)"
	echo "    -or | --overlap_read INT:         Minimum length of overlap for creating edges between two read in the read graph. (default: 5)"
	echo "    -wr | --weight_read FLOAT:        Minimum weights of edges in the read graph. (default:0.8)"
	echo "    -m  | --mcl_inflaction FLOAT:     Inflaction of MCL algorithm. (default:2)"
	echo "    -l  | --lar_cluster INT:          A threshold for seperating clusters into two groups based on sizes of clusters. (default:50)"
	echo "    -oc | --overlap_cluster INT:      A parameter related to the minimum overlap between consensus sequences. (default:10) "
	echo "    -d  | --depth INT:                Depth limitation for consensus sequences generated from clusters. (default:5) "
	echo "    -wc | --weight_cluster FLOAT:     Minimum weights between clusters in the hierarchical clustering. (default: 0.8)"
	echo "    -a  | --abundance FLOAT:          A threshold for filtering low-abundance haplotypes. (default: 0.005)"
	echo "    -h  | --help :                    Print help message."
	echo ""
	echo "    For further details of above arguments, please refer to https://github.com/dhcai21/RVHaplo   "
	echo ""
	exit 1
}

if [[ "$1" == "" ]];then
	help_info
	exit 1
fi

while [[ "$1" != "" ]]; do
	case "$1" in
		-h | --help ) ## print help message
		help_info
		exit 1
		;;
		-i | --input ) ### input sam file
		case "$2" in 
		"" )
			echo "Error: no sam file as input"
			exit 1
			;;
		*)
			file_sam="$2"
			if [[ "${file_sam:0:1}" == "-" ]]
			then
				echo "Error: no sam file as input"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-r | --ref_genome) ### input reference genome
		case "$2" in 
		"")
			echo "Error: no fasta file as input"
			exit 1
			;;
		*)
			file_ref="$2"
			if [[ ""${file_ref:0:1}"" == "-" ]]
			then
				echo "Error: no fasta file as input"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-o | --out )  ### output path
		case "$2" in 
		"" )
			echo "Error: no output path"
			exit 1
			;;
		*)
			file_path="$2"
			if [[ "${file_sam:0:1}" == "-" ]]
			then
				echo "Error: no output path"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-p | --prefix )  ### prefix
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			prefix="$2"
			shift 2
			;;
		esac
		;;
		-e | --error_rate )  ### error_rate
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			error_rate="$2"
			shift 2
			;;
		esac
		;;
		-s | --signi_level )  ### significance_level
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			signi_level="$2"
			shift 2
			;;
		esac
		;;
		-c | --cond_pro )  ### conditional_probability
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			cond_pro="$2"
			shift 2
			;;
		esac
		;;
		-f | --fre_snv )  ### determine the set of to-be-verified SNV sites
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			fre_snv="$2"
			shift 2
			;;
		esac
		;;
		-n1 | --num_read_1 )  ### number of reads for p(ai|aj)
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			num_read_1="$2"
			shift 2
			;;
		esac
		;;
		-n2 | --num_read_2 )  ### number of reads for p(ai|aj1,aj2,...,ajp)
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			num_read_2="$2"
			shift 2
			;;
		esac
		;;
		-g | --gap )  ### Minimum distance between SNV sites
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			gap="$2"
			shift 2
			;;
		esac
		;;
		-s | --smallest_snv )  ### Minimum number of SNV sites for haplotype reconstruction
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			smallest_snv="$2"
			shift 2
			;;
		esac
		;;
		-or | --ovlap_read )  ### overlap_read
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			ovlap_read="$2"
			shift 2
			;;
		esac
		;;
		-wr | --weight_read )  ### weight_read
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			weight_read="$2"
			shift 2
			;;
		esac
		;;
		-m | --mcl_inflaction )  ### inflaction of MCL
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			mcl_inflaction="$2"
			shift 2
			;;
		esac
		;;
		-oc | --ovlap_cluster )  ### overlap_cluster
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			ovlap_cluster="$2"
			shift 2
			;;
		esac
		;;
		-d | --depth )  ### Minimum depth of consensus sequence
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			ovlap_cluster="$2"
			shift 2
			;;
		esac
		;;
		-wc | --weight_cluster )  ### weight_cluster
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			weight_cluster="$2"
			shift 2
			;;
		esac
		;;
		-d | --depth )  ### depth limitation
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			depth="$2"
			shift 2
			;;
		esac
		;;
		-l | --lar_cluster )  ### large cluster size
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			lar_cluster="$2"
			shift 2
			;;
		esac
		;;
		-a | --abundance )  ### smallest abundance
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			abundance="$2"
			shift 2
			;;
		esac
		;;
		*)
			echo "Error: unknow parameter $1"
			exit 1
	esac
done

if [[ "$file_sam" == "" ]];then
	echo "Error: no sam file input"
	exit 1
fi

if [[ "$file_ref" == "" ]];then
	echo "Error: no reference genome input"
	exit 1
fi

if [[ ${file_path:0-1:1} == "/" ]];then
	path_len=`expr ${#file_path}-1`
	file_prefix=$file_path$prefix
	file_path=${file_path:0:path_len}
else
	file_prefix=$file_path"/"$prefix
fi

##########  count nucleotide occurrence  ##########
echo "count nucleotide occurrence"
file_len=`expr ${#file_sam}-4`
unique_sam=${file_sam:0:$file_len}"_unique.sam"
samtools view -h -F 0x900 $file_sam > $unique_sam
file_bam=${file_sam:0:$file_len}".bam"
samtools view -b -S $unique_sam > $file_bam
file_bam_sorted=${file_sam:0:$file_len}"_sorted.bam"
samtools sort $file_bam -o $file_bam_sorted
samtools index $file_bam_sorted
file_acgt=$file_prefix"_acgt.txt"
pysamstats -f $file_ref --type variation_strand $file_bam_sorted > $file_acgt


########## two binomial tests  ##########
echo "SNV detection"
file_snv=$file_prefix"_snv.txt"
python ./src/two_binomial.py $error_rate $signi_level $file_acgt $file_snv

## judge number of detected SNV sites
size="$(wc -l <"$file_snv")"
if [[ $size != "0" ]];then
	exit 1
fi

## maximum conditional probability and apply MCL to read graph 
echo "Graph clustering"
python ./src/read_graph_mcl.py $file_bam_sorted $file_snv $cond_pro $smallest_snv $num_read_1 $num_read_2 $gap \
	$weight_read $ovlap_read $mcl_inflation $file_prefix $fre_snv

## judge number of detected SNV sites
size="$(wc -l <"$file_snv")"
if [[ $size != "0" ]];then
	exit 1
fi

## hierarchical clustering
echo "hierarchical clustering"
python ./src/hierarchical_cluster.py $file_prefix"_graph.pickle" $lar_cluster $depth \
	$ovlap_cluster $weight_cluster $abundance $file_prefix

## reconstruct haplotypes
rm -rf $file_path"/clusters"
mkdir $file_path"/clusters"

echo "haplotypes reconstruction"

python ./src/out_haplotypes.py $file_prefix"_final.pickle" $file_bam $file_path $file_acgt $file_ref \
	$file_prefix"_haplotypes.fasta"

echo "complete reconstructing haplotypes"

exit 1
