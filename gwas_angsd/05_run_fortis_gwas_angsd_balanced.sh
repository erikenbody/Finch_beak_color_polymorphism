#! /bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 6
#SBATCH -t 1-00:00:00
#SBATCH -J angsd_bal_fortis
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

if [ -d "Results_gwas_BALENCED" ]; then echo "Results file exists" ; else mkdir Results_gwas_BALENCED; fi

POP1=fortis_yellow
POP2=fortis_pink

##do once
# shuf -n 130 ${POP2}_bams.txt > ${POP2}_bams_SUBSET.txt
#
# cat ${POP1}_bams.txt ${POP2}_bams_SUBSET.txt > ${POP1}_${POP2}_bams_BALENCED.txt
#
# #create pheno:
# printf '1\n%.0s' {1..130} > fortis_bin_pheno_BALENCED.txt
# printf '0\n%.0s' {1..130} >> fortis_bin_pheno_BALENCED.txt

BAMLIST1=${POP1}_${POP2}_bams_BALENCED.txt

$ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $REFGENOME -r $INTERVAL \
              -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref \
              -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
              -minMapQ 20 -minQ 20 -doCounts 1 \
              -domajorminor 1 -domaf 1 \
              -GL 1 -P 6 -doGlf 2 -SNP_pval 1e-6 -minMaf 0.05 \
              -dumpCounts 1

#have to run pc, create phenotype file
python $PCANGSD/pcangsd/pcangsd.py -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.beagle.gz -o Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.pcangsd -threads 6
Rscript ~/bc/gwas_angsd/get_PCs.R Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.pcangsd.cov

#then will need to make a covariate file for doass

# $ANGSD_PATH/angsd -doMaf 4 -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 6 -cov Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.res
# $ANGSD_PATH/angsd -Pvalue 1 -doMaf 4 -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 6 -cov Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.res.p
# $ANGSD_PATH/angsd -Pvalue 1 -doMaf 4 -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 4 -cov Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.res.p.4
# $ANGSD_PATH/angsd -Pvalue 1 -doMaf 4 -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 2 -model recessive -cov Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.res.p.2.rec
# $ANGSD_PATH/angsd -Pvalue 1 -doMaf 4 -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 6 -model recessive -cov Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.res.p.6.rec
# $ANGSD_PATH/angsd -doMaf 4 -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno.txt -doAsso 6 -model recessive -cov Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.ref.PC1_PC2.txt -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}.res.lrt.6.rec

#score test (2) and EM model (4). Not sure the value of dosage (it basically calls genotypes). LRT and Pvalues provide different results, but LRT seems more consistent across cases
#recessive model works

$ANGSD_PATH/angsd -doMaf 4 -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno_BALENCED.txt -doAsso 2 -model recessive -cov Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.PC1_PC2.txt -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.lrt.2.rec
#$ANGSD_PATH/angsd -doMaf 4 -beagle Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.beagle.gz -fai $FAIFILE  -yBin fortis_bin_pheno_BALENCED.txt -doAsso 4 -model recessive -cov Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.PC1_PC2.txt -out Results_gwas_BALENCED/${POP1}_${POP2}_${INTERVAL}_BALENCED.ref.lrt.4.rec
