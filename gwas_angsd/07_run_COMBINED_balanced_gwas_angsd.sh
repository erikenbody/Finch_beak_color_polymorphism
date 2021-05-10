#! /bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 8 #was 6 cores for all but 1 and 2
#SBATCH -t 7-00:00:00
#SBATCH -J angsd_COMBINED
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e angsd_%A_%a.err            # File to which STDERR will be written
#SBATCH -o angsd_%A_%a.out

module load samtools/1.10

ANGSD_PATH=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/angsd_14April/angsd
PCANGSD=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/pcangsd_15Apr/
FAIFILE=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta.fai
REFGENOME=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta
INTERVAL_FILE=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/Cam_scaffolds_corrected.txt

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_BCO2/angsd/revision_2021
cd $WORK_D

INTERVAL=`cat $INTERVAL_FILE | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
echo $INTERVAL
#24 and chr24 is BC02

if [ -d "Results_gwas" ]; then echo "Results file exists" ; else mkdir Results_gwas; fi

POP1=combined_yellow
POP2=combined_pink

#cat scandens_yellow_scandens_pink_bams_BALENCED.txt fortis_yellow_fortis_pink_bams_BALENCED.txt > combined_yellow_combined_pink_BALENCED_bams.txt
#cat scandens_bin_pheno_BALENCED.txt fortis_bin_pheno_BALENCED.txt > combined_bin_pheno_BALENCED.txt

BAMLIST1=${POP1}_${POP2}_BALENCED_bams.txt
#p6 for most
$ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $REFGENOME -r $INTERVAL \
              -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref \
              -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
              -minMapQ 20 -minQ 20 -doCounts 1 \
              -domajorminor 1 -domaf 1 \
              -GL 1 -P 8 -doGlf 2 -SNP_pval 1e-6 -minMaf 0.05 \
              -dumpCounts 1

#have to run pc, create phenotype file
python $PCANGSD/pcangsd/pcangsd.py -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -o Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.pcangsd -threads 8
Rscript ~/bc/gwas_angsd/get_PCs.R Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.pcangsd.cov

#then will need to make a covariate file for doass

# $ANGSD_PATH/angsd -doMaf 4 -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 6 -cov Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.res
# $ANGSD_PATH/angsd -Pvalue 1 -doMaf 4 -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 6 -cov Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.res.p
# $ANGSD_PATH/angsd -Pvalue 1 -doMaf 4 -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 4 -cov Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.res.p.4
# $ANGSD_PATH/angsd -Pvalue 1 -doMaf 4 -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 2 -model recessive -cov Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.res.p.2.rec
# $ANGSD_PATH/angsd -Pvalue 1 -doMaf 4 -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 6 -model recessive -cov Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.res.p.6.rec
# $ANGSD_PATH/angsd -doMaf 4 -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 6 -model recessive -cov Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.res.lrt.6.rec

#score test (2) and EM model (4). Not sure the value of dosage (it basically calls genotypes). LRT and Pvalues provide different results, but LRT seems more consistent across cases
#recessive model works

$ANGSD_PATH/angsd -doMaf 4 -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin combined_bin_pheno_BALENCED.txt -doAsso 2 -model recessive -cov Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.res.lrt.2.rec
#$ANGSD_PATH/angsd -doMaf 4 -beagle Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 4 -model recessive -cov Results_gwas/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas/${POP1}_${POP2}_${INTERVAL}.res.lrt.4.rec
