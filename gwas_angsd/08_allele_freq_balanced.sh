#! /bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 3
#SBATCH -t 1-00:00:00
#SBATCH -J angsd
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e angsd_%A_%a.err            # File to which STDERR will be written
#SBATCH -o angsd_%A_%a.out

module load samtools/1.10

ANGSD_PATH=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/angsd_14April/angsd
PCANGSD=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/pcangsd_15Apr/
FAIFILE=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta.fai
REFGENOME=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta
ANCESTRAL=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_BCO2/angsd/revision_2021/Results_fortis_scandens/ancestral_chr24.fa.gz
INTERVAL_FILE=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/Cam_scaffolds_corrected.txt

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_BCO2/angsd/revision_2021
cd $WORK_D

INTERVAL=chr24
echo $INTERVAL
#24 and chr24 is BC02

if [ -d "Results_af" ]; then echo "Results file exists" ; else mkdir Results_af; fi
POP1=fortis_yellow
POP2=fortis_pink_bams_SUBSET
POP3=scandens_yellow
POP4=scandens_pink_bams_SUBSET


SITES=top_snps_fortis_scandens_balanced.txt
$ANGSD_PATH/angsd sites index $SITES

##polarize alleles with the ancestral allele and only take high confidence variable sites from the right hand of chr24
BAMLIST1=${POP1}_bams.txt
$ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $ANCESTRAL -r $INTERVAL -sites $SITES \
              -out Results_af/${POP1}_${INTERVAL}_BALANCED.ref \
              -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
              -minMapQ 20 -minQ 20 -doCounts 1 \
              -domajorminor 4 -domaf 2 \
              -GL 1 -P 5 -SNP_pval 1e-6

BAMLIST1=${POP2}.txt

$ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $ANCESTRAL -r $INTERVAL -sites $SITES \
              -out Results_af/${POP2}_${INTERVAL}_BALANCED.ref \
              -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
              -minMapQ 20 -minQ 20 -doCounts 1 \
              -domajorminor 4 -domaf 2 \
              -GL 1 -P 5 -SNP_pval 1e-6

BAMLIST1=${POP3}_bams.txt

$ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $ANCESTRAL -r $INTERVAL -sites $SITES \
              -out Results_af/${POP3}_${INTERVAL}_BALANCED.ref \
              -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
              -minMapQ 20 -minQ 20 -doCounts 1 \
              -domajorminor 4 -domaf 2 \
              -GL 1 -P 5 -SNP_pval 1e-6

BAMLIST1=${POP4}.txt

$ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $ANCESTRAL -r $INTERVAL -sites $SITES \
              -out Results_af/${POP4}_${INTERVAL}_BALANCED.ref \
              -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
              -minMapQ 20 -minQ 20 -doCounts 1 \
              -domajorminor 4 -domaf 2 \
              -GL 1 -P 5 -SNP_pval 1e-6
