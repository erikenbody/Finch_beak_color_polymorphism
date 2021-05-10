#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -p core -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J allele_sharing
#SBATCH -e allele_sharing_%A_%a.err            # File to which STDERR will be written
#SBATCH -o allele_sharing_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.10

WORK_D=/home/eenbody/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_BCO2/allele_sharing2
cd $WORK_D

VCF=/home/eenbody/snic2020-2-19/private/darwins_finches/variants/FILTERED/293.SNP.vcf.gz

grep -v "Mixed" sample_pop_ID.txt  | grep -v "vcf" | grep -v "crass" | grep -v "olivacea" | grep -v "fusca" | grep -v "BB" | grep -v "noctis" | grep -v "x" | grep -v "bicolor" | uniq > sample_pop_ID_uniq.txt
sort -u -k1,1 sample_pop_ID_uniq.txt | cut -f 1 > species.txt

##only ground and tree with 90% samples have data
bcftools view -S vcf_samps_include.txt -e 'AC==0 || AC==AN || F_MISSING > 0.1 || ALT="*"' -m2 -M2 -O b -o ground_tree_missing0.1_bial.bcf $VCF
bcftools index ground_tree_missing0.1_bial.bcf

#need to run with AC filters again after sub sampling the file
bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.1 || ALT="*"' -m2 -M2 -O b -o ground_tree_missing0.1_bial_fil.bcf ground_tree_missing0.1_bial.bcf
bcftools index ground_tree_missing0.1_bial_fil.bcf

BCF=ground_tree_missing0.1_bial_fil.bcf

#calculate allele frequency per species
while read species
do
  grep $species sample_pop_ID_uniq.txt | cut -f 2 > ${species}.list
  echo $species
  bcftools view -S ${species}.list $BCF | bcftools query -f '%CHROM %POS %AN %AC{0}\n' > ${species}.af
  cut -f 4 -d  ' ' ${species}.af | awk '$1!=0 {$1=1} {print}' > ${species}.bin.af
done < species.txt

paste *bin.af > all_species.af

awk '{ for(i=1; i<=NF;i++) j+=$i; print j; j=0 }' all_species.af > sum_alt_all_species.txt

while read species
do
  sed -n 17p ${species}.af
done < species.txt
