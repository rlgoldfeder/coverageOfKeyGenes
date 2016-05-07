#!/bin/sh

################################################
# author: Rachel Goldfeder
# Calculate Coverage (or lack thereof) For Key Genes
# Requires bedtools
################################################




usage ()
{
  echo 'Usage : coverageOfKeyGenes.sh -gatk <gatk_path> -bq <bq> -mq <mq>' 
  echo '                              -cov <cov> -g <geneList.bed> -bam <bam>'
  echo '                              -prefix <bam_prefix> -r <ref.fa> -mem <gatk_gigs_mem>'
  exit
}

if [ "$#" -ne 8 ]
then
  usage
fi

while [ "$1" != "" ]; do
case $1 in
        -gatk )           shift
                       gatk_path=$1
                       ;;
        -bq )           shift
                       bq=${1:-20}
                       ;;
        -mq )           shift
                       mq=${1:-1}
		       ;;
        -cov )           shift
                       cov=${1:-20}
                       ;;
       	 -g )           shift
                       gene_list_path=$1
                       ;;
        -bam )           shift
                       bam_path=$1
                       ;;
        -prefix )           shift
                       bam_prefix=$1
                       ;;
        -r )           shift
                       ref_genome=$1
                       ;;
	-mem )           shift
                       mem=${1:-16}
                       ;;
        * )            QUERY=$1
    esac
    shift
done


################################################
# 1. Run GATK Depth of Coverage
################################################

mkdir coverage

java -Xmx"$mem"G -jar "$gatk_path" \
	  -T DepthOfCoverage \
	  -mbq "$bq" \
	  -mmq "$mq" \
	  -L "$gene_list_path" \
	  -I "$bam_path" \
	  --summaryCoverageThreshold 0 \
	  -R "$ref_genome" \
	  -o coverage/"$bam_prefix"."$bq"bq."$mq"mq.txt 

################################################
# 2. Turn GATK output text file into a bed file
################################################

awk -v OFS="\t" -v c="$cov" 'NR>1{if($2 >= c){split($1, a, ":"); print a[1],a[2]-1,a[2]}}' coverage/"$bam_prefix"."$bq"bq."$mq"mq.txt  \
	  >  coverage/"$bam_prefix"."$bq"bq."$mq"mq.bed


################################################
# 3. Intersect output bed with orig gene list
################################################
words=$(wc -l < coverage/"$bam_prefix"."$bq"bq."$mq"mq.bed)

if [ $words -eq 0 ];then
	 # echo "None of the bases in any of your genes pass your thresholds. Make sure your gene coordinates are correct or try more lenient thresholds."
	  #  exit 1
	   
	  awk -v OFS="\t" '{arr[$4]+=0} END {for (i in arr) {print i, arr[i]}}' $gene_list_path > coverage/"$bam_prefix"."$bq"bq."$mq"mq.genes.sum.bed  

else
	 bedtools intersect  \
    	-a $gene_list_path \
	-b coverage/"$bam_prefix"."$bq"bq."$mq"mq.bed \
        -wao \
	 > coverage/"$bam_prefix"."$bq"bq."$mq"mq.genes.bed


################################################
# 4. Sum By Gene
################################################

awk -v OFS="\t" '{arr[$4]+=$8} END {for (i in arr) {print i, arr[i]}}' coverage/"$bam_prefix"."$bq"bq."$mq"mq.genes.bed > coverage/"$bam_prefix"."$bq"bq."$mq"mq.genes.sum.bed
fi

sort coverage/"$bam_prefix"."$bq"bq."$mq"mq.genes.sum.bed >  coverage/"$bam_prefix"."$bq"bq."$mq"mq.genes.sum.sorted.bed
awk -v OFS="\t" '{arr[$4]+=$3-$2} END {for (i in arr) {print i,arr[i]}}' $gene_list_path  > coverage/totByGene.bed
sort coverage/totByGene.bed > coverage/sorted.totByGene.bed
		  
join -1 1 -2 1 -t $'\t' coverage/"$bam_prefix"."$bq"bq."$mq"mq.genes.sum.sorted.bed coverage/sorted.totByGene.bed > coverage/"$bam_prefix"."$bq"bq."$mq"mq.summary.txt
		    
		      
awk -v OFS="\t" 'BEGIN {print "Gene\tnumMissing\tfracMissing"} {tot_fail = $3-$2; frac_fail = tot_fail / $3; print $1,tot_fail,frac_fail }' coverage/"$bam_prefix"."$bq"bq."$mq"mq.summary.txt > coverage/"$bam_prefix"."$bq"bq."$mq"mq."$cov".cov.coverage_metrics.txt

 
 
################################################
# 5. Plot results
################################################

Rscript --vanilla plot_coverage.R coverage/"$bam_prefix"."$bq"bq."$mq"mq."$cov".cov.coverage_metrics.txt $cov $bq $mq




