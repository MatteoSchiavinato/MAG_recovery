#!/usr/bin/env sh

PROJECT="M00001_PacBio_chicken_gut_metagenomics"
WD="/gpfs/data/fs71579/schmat90/CF/projects/meta/${PROJECT}"

echo """\
#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -J MP-MAG
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user matteo.schiavinato@boku.ac.at
#SBATCH --account p71579
#SBATCH --qos mem_0384
#SBATCH --partition mem_0384

module purge
module load gcc/10.2.0-gcc-9.1.0-2aa5hfe

cd ${WD}/scripts/wf-MAG_recovery

nextflow \
run \
main.nf \
-resume \
-work-dir ${WD}/scripts/wf-MAG_recovery/work \
-with-report ${WD}/scripts/wf-MAG_recovery/cmd.sbatch.report.html \
-with-timeline ${WD}/scripts/wf-MAG_recovery/cmd.sbatch.timeline.html \
-with-dag ${WD}/scripts/wf-MAG_recovery/cmd.sbatch.dag.png \
--output_dir ${WD} \
--threads 48 \
--pplacer_threads 8 \
--min_cluster_tot_len 100000 \
--contigs_fasta ${WD}/assembly/canu/ALL.ccs.contigs.fasta \
--ccs_reads ${WD}/raw_data/demux_ccs_fasta/ALL.fasta \
--sample_id ALL \
--min_perc_aa 10 \
--min_align_frac 0.65 \

""" \
> cmd.sbatch

sbatch --output cmd.sbatch.out cmd.sbatch
