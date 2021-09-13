#! /bin/bash
#SBATCH -A snic2021-5-8
#SBATCH -p node #node for 1-5, 8 for others
#SBATCH -t 6-00:00:00
#SBATCH -J angsd
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e angsd_%A_%a.err            # File to which STDERR will be written
#SBATCH -o angsd_%A_%a.out

#following getting negative Fst values and reading this:
#https://github.com/ANGSD/angsd/issues/274

module load samtools/1.10

ANGSD_PATH=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/angsd_14April/angsd
WORK_D=//crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_BCO2/angsd/revision_2021/fst_species
REFGENOME=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta

INTERVAL_FILE=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/Cam_scaffolds_corrected.txt

INTERVAL=`cat $INTERVAL_FILE | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
echo $INTERVAL
#24 and chr24 is BC02

cd $WORK_D

if [ -d "Results_fortis_scandens" ]; then echo "Results file exists" ; else mkdir Results_fortis_scandens; fi

POP1=fortis_GG_early
POP2=scandens_GG_early
POP3=fortis_GG_late
POP4=scandens_GG_late

BAMLIST1=${POP1}_bams.txt

$ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $REFGENOME -r $INTERVAL \
              -out Results_fortis_scandens/${POP1}_${INTERVAL}.ref \
              -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
              -minMapQ 20 -minQ 20 -doCounts 1 \
              -GL 1 -P 8 \
              -doSaf 1

# BAMLIST1=${POP2}_bams.txt
#
# $ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $REFGENOME -r $INTERVAL \
#               -out Results_fortis_scandens/${POP2}_${INTERVAL}.ref \
#               -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
#               -minMapQ 20 -minQ 20 -doCounts 1 \
#               -GL 1 -P 8 \
#               -doSaf 1
#
# BAMLIST1=${POP3}_bams.txt
#
# $ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $REFGENOME -r $INTERVAL \
#               -out Results_fortis_scandens/${POP3}_${INTERVAL}.ref \
#               -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
#               -minMapQ 20 -minQ 20 -doCounts 1 \
#               -GL 1 -P 8 \
#               -doSaf 1
#
# BAMLIST1=${POP4}_bams.txt
#
# $ANGSD_PATH/angsd -b $BAMLIST1 -ref $REFGENOME -anc $REFGENOME -r $INTERVAL \
#               -out Results_fortis_scandens/${POP4}_${INTERVAL}.ref \
#               -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
#               -minMapQ 20 -minQ 20 -doCounts 1 \
#               -GL 1 -P 8 \
#               -doSaf 1

#re running with the proper folding commands
#
$ANGSD_PATH/misc/realSFS -P 20 Results_fortis_scandens/${POP1}_${INTERVAL}.ref.saf.idx Results_fortis_scandens/${POP2}_${INTERVAL}.ref.saf.idx > Results_fortis_scandens/${POP1}_${POP2}_${INTERVAL}_unfolded.sfs
$ANGSD_PATH/misc/realSFS -P 20 Results_fortis_scandens/${POP3}_${INTERVAL}.ref.saf.idx Results_fortis_scandens/${POP4}_${INTERVAL}.ref.saf.idx > Results_fortis_scandens/${POP3}_${POP4}_${INTERVAL}_unfolded.sfs
$ANGSD_PATH/misc/realSFS -P 20 Results_fortis_scandens/${POP1}_${INTERVAL}.ref.saf.idx Results_fortis_scandens/${POP4}_${INTERVAL}.ref.saf.idx > Results_fortis_scandens/${POP1}_${POP4}_${INTERVAL}_unfolded.sfs

$ANGSD_PATH/misc/realSFS fst index Results_fortis_scandens/${POP1}_${INTERVAL}.ref.saf.idx Results_fortis_scandens/${POP2}_${INTERVAL}.ref.saf.idx -sfs Results_fortis_scandens/${POP1}_${POP2}_${INTERVAL}_unfolded.sfs -fstout Results_fortis_scandens/${POP1}_${POP2}_${INTERVAL}_unfolded
$ANGSD_PATH/misc/realSFS fst index Results_fortis_scandens/${POP3}_${INTERVAL}.ref.saf.idx Results_fortis_scandens/${POP4}_${INTERVAL}.ref.saf.idx -sfs Results_fortis_scandens/${POP3}_${POP4}_${INTERVAL}_unfolded.sfs -fstout Results_fortis_scandens/${POP3}_${POP4}_${INTERVAL}_unfolded
$ANGSD_PATH/misc/realSFS fst index Results_fortis_scandens/${POP1}_${INTERVAL}.ref.saf.idx Results_fortis_scandens/${POP4}_${INTERVAL}.ref.saf.idx -sfs Results_fortis_scandens/${POP1}_${POP4}_${INTERVAL}_unfolded.sfs -fstout Results_fortis_scandens/${POP1}_${POP4}_${INTERVAL}_unfolded

$ANGSD_PATH/misc/realSFS fst stats2 Results_fortis_scandens/${POP1}_${POP2}_${INTERVAL}_unfolded.fst.idx -win 10000 -step 10000 -whichFST 0 > Results_fortis_scandens/${POP1}_${POP2}__${INTERVAL}_unfolded_10kb-Win.fst.txt
$ANGSD_PATH/misc/realSFS fst stats2 Results_fortis_scandens/${POP3}_${POP4}_${INTERVAL}_unfolded.fst.idx -win 10000 -step 10000 -whichFST 0 > Results_fortis_scandens/${POP3}_${POP4}__${INTERVAL}_unfolded_10kb-Win.fst.txt
$ANGSD_PATH/misc/realSFS fst stats2 Results_fortis_scandens/${POP1}_${POP4}_${INTERVAL}_unfolded.fst.idx -win 10000 -step 10000 -whichFST 0 > Results_fortis_scandens/${POP1}_${POP4}__${INTERVAL}_unfolded_10kb-Win.fst.txt

INTERVAL=chr24
$ANGSD_PATH/misc/realSFS fst print Results_fortis_scandens/${POP1}_${POP2}_${INTERVAL}_unfolded.fst.idx -r chr24:6156878-6176878 > ${POP1}_${POP2}_${INTERVAL}_5kb_region.fst
$ANGSD_PATH/misc/realSFS fst print Results_fortis_scandens/${POP3}_${POP4}_${INTERVAL}_unfolded.fst.idx -r chr24:6156878-6176878 > ${POP3}_${POP4}_${INTERVAL}_5kb_region.fst
