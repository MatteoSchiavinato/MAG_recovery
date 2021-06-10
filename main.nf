#!/usr/bin/env nextflow

// help section
if (params.help) {

  println """

  Matteo schiavinato
  Metagenomic contig clustering and discovery of MAGs
  March 2021-June 2021

  ### misc ###

  --threads               Number of parallel threads
  --sample_id             Name of the processed sample

  ### input ###

  --contigs_fasta         FASTA file containing assembled contigs
  --ccs_reads             FASTA file containing original ccs reads

  ### contig clustering into bins ###

  --min_contig_length	    Min. length of a contig to retain in bin (1500)
  --perc_good_contigs 	  % of contigs in bin that share similar ACTG-freq (95)
  --min_ctg_coverage	    Min. coverage of contigs (1)
  --min_cluster_tot_len	  Min. total length of a bin for it to be retained (100000)

  ### bin lineage and filtering ###

  --pplacer_threads	      Threads to be used by the pplacer algorithm, used in "checkm tree" (16)
  --max_contamination     Maximum % of contamination to retain a bin, from 0 to 100 (5.0)
  --min_completeness      Minimum % of completeness to retain a bin, from 0 to 100 (10.0)
  --max_heterogeneity     Maximum % of heterogeneity to retain a bin, from 0 to 100 (25.0)

  ### bin refinement ###

  --minimap_kmer_size     Kmer size to use for minimap2 mapping, to extract coverage (20)
  --min_align_len         Minimum fraction of ccs read that has to map to consider it for coverage (0.75)
  --cov_max_edit_dist     Maximum edit distance (frac) for a ccs read to consider it for coverage (0.10)
  --gc_perc               Minimum percentile to consider scaffold not outlier in terms of GC content (95)
  --td_perc               Minimum percentile to consider scaffold not outlier in terms of tetranucleotide freq (95)
  --cov_perc              Minimum percentile to consider scaffold not outlier in terms of coverage (95)

  ### bin taxa annotation ###

  --min_perc_aa           Minimum % of aminoacids that a MAG must contain in the MSA to be retained (10)
  --min_align_frac        Minimum fraction of MAG aligned to consider hit to a closest genome (0.65)

  ### output ###

  --output_dir            Base directory where all the sub-directories with the outputs
                          will be placed.
                          Ideally, this is the directory containing /raw_data and /scripts

  ### more about metabat2 usage:
  ### https://bitbucket.org/berkeleylab/metabat/wiki/Home

  ### more about CheckM usage:
  ### http://ecogenomics.github.io/CheckM/manual/checkm_manual.pdf

  ### more about RefineM usage:
  ### https://github.com/dparks1134/RefineM

  ### more about gtdb-tk usage:
  ### https://ecogenomics.github.io/GTDBTk/index.html

  """
  exit 0
}

// -----------------------------------------------------------------------------

// binning and generation of MAGs


// checkm Manual:
// http://ecogenomics.github.io/CheckM/manual/checkm_manual.pdf

// binning assemblies

process metabat2 {

  executor = 'local'
  cpus = params.threads

  publishDir "${params.output_dir}/metabat2", mode: "copy"

  output:
    path("${params.sample_id}_out", hidden:false, type:"dir") into Metabat2_out

  script:
    """
    ${METABAT2} \
    --inFile ${params.contigs_fasta} \
    --outFile ${params.sample_id}_out/${params.sample_id}.bin \
    --minContig ${params.min_contig_length} \
    --maxP ${params.perc_good_contigs} \
    --minS ${params.min_edge_score} \
    --maxEdges ${params.max_num_edges} \
    --minCV ${params.min_ctg_coverage} \
    --minClsSize ${params.min_cluster_tot_len} \
    --numThreads ${params.threads}
    """
}

// -----------------------------------------------------------------------------
// WARNING:
// before running any checkm command you have to run: checkm data setRoot db/
// where db/ is the database you downloaded
// this is documented upon installation of checkm
// but I thought of reporting it here too, just so that whoever is reading knows
// -----------------------------------------------------------------------------

// to use checkm properly, consult their wiki:
// https://github.com/Ecogenomics/CheckM/wiki

Metabat2_out.into{Bins_for_filtering; Bins_for_lineage_wf}



// running checkM with lineage-specific parameters
// filtering for completeness and contamination
// this is based on a set of marker genes
// in paper: 70% completeness, 5% contamination
// running lineage-wf

process checkm_lineage_wf {

  executor = 'local'
  cpus = params.threads

  publishDir "${params.output_dir}/checkm", mode: "copy"

  input:
    path bins from Bins_for_lineage_wf

  output:
    path "${params.sample_id}_raw_lineage", type: "dir" into Checkm_lineage_out
    file "checkm_lineage_wf_results.txt" into lineage_file
    file "${params.sample_id}_raw_lineage/lineage.ms" into marker_file

  script:
    """
    if [ ! -d ${params.sample_id}.TMP ] ; then mkdir ${params.sample_id}.TMP ; fi
    ${CHECKM} \
    lineage_wf \
    --ali \
    --nt \
    --e_value ${params.e_value} \
    --length ${params.overlap_length_frac} \
    --file checkm_lineage_wf_results.txt \
    --threads ${params.threads} \
    --extension fa \
    --pplacer_threads ${params.pplacer_threads} \
    ${bins} \
    ${params.sample_id}_raw_lineage \

    """
}


// filter bins based on contamination and completeness

process filter_bins {

  executor = 'local'
  cpus = params.threads

  publishDir "${params.output_dir}/checkm", mode: "copy"

  input:
    file lineage_file
    path bins from Bins_for_filtering

  output:
    file "checkm_lineage_wf_results.tsv"
    file "filter_bins.log"
    path "${params.sample_id}_filtered_lineage"
    path "${params.sample_id}_filtered_lineage/bins" into Filtered_bins

  script:
    """
    { echo -e "Bin_Id\tMarker_lineage\t#_genomes\t#_markers\t#_marker_sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain_heterogeneity"; \
    tail -n+4 ${lineage_file} | \
    sed 's/ (/_(/' | \
    tr -s " " "\t" | \
    cut -f 2- | \
    awk 'NF != 1'; } \
    > checkm_lineage_wf_results.tsv &&
    ${PYTHON3} \
    ${params.source_dir}/filter-bins.py \
    --max-cont ${params.max_contamination} \
    --min-comp ${params.min_completeness} \
    --max-het ${params.max_heterogeneity} \
    --bins-directory ${bins} \
    --extension fa \
    --lineage-file checkm_lineage_wf_results.tsv \
    --output-directory ${params.sample_id}_filtered_lineage \
    2> filter_bins.log \

    """
}

Filtered_bins.into{ Bins_filt_merge;
                    Bins_filt_unique;
                    Bins_filt_refinem_stats;
                    Bins_filt_filter_bins; }



// ensuring no contig is assigned to multiple bins

process checkm_unique {

  executor = 'local'
  cpus = 1

  publishDir "${params.output_dir}/statistics", mode: "copy"

  input:
  path bins from Bins_filt_unique

  output:
  file "checkm_unique_out.txt"

  script:
  """
  ${CHECKM} \
  unique \
  -x fa \
  ${bins} \
  &> checkm_unique_out.txt
  """
}

// map ccs reads to get coverage profile
// using minimap2

process ccs_mapping {

  executor = 'local'
  cpus = 1

  publishDir "${params.output_dir}/ccs_mapping", mode: "copy"

  output:
    file "${params.sample_id}.bam" into ccs_bam

  script:
    """
    ${MINIMAP2} \
    -H \
    -k ${params.minimap_kmer_size} \
    -d contigs.mm2_index \
    -a \
    -t ${params.threads} \
    -x map-pb \
    ${params.contigs_fasta} \
    ${params.ccs_reads} \
    > ${params.sample_id}.bam \

    """
}


// filter and sort the bam file
// using samtools

process ccs_map_filter_sort {

  executor = 'local'
  cpus = 1

  publishDir "${params.output_dir}/ccs_mapping", mode: "copy"

  input:
    file ccs_bam

  output:
    file "${params.sample_id}.f.bam" into ccs_filt_bam
    file "${params.sample_id}.fs.bam" into ccs_filt_sort_bam

  script:
    """
    ${SAMTOOLS} \
    view \
    -b -h -@ ${params.threads} \
    ${ccs_bam} \
    -F 0x0100 -F 0x4 \
    > ${params.sample_id}.f.bam &&
    ${SAMTOOLS} \
    sort \
    -@ ${params.threads} \
    -T ${params.sample_id} \
    ${params.sample_id}.f.bam \
    > ${params.sample_id}.fs.bam \

    """
}


// filter bam file
// using samtools

process index_bam {

  input:
    file ccs_filt_sort_bam

  output:
    file "${ccs_filt_sort_bam}.bai" into ccs_filt_sort_bam_index

  script:
    """
    ${SAMTOOLS} \
    index \
    ${ccs_filt_sort_bam} \

    """
}


// extract unmapped reads to make statistics
// with samtools

process recover_unmapped_reads {

  executor = 'local'
  cpus = 1

  publishDir "${params.output_dir}/ccs_mapping", mode: "copy"

  input:
    file ccs_filt_bam

  output:
    file "${params.sample_id}.unmapped.fasta" into unmapped_reads_fasta

  script:
    """
   ${SAMTOOLS} fasta \
   -@ ${params.threads} \
   ${ccs_filt_bam} \
   > ${params.sample_id}.unmapped.fasta \

    """
}


// detect scaffold statistics
// with refinem scaffold_stats
// and the bam file from ccs to get coverage from

process refinem_scaffold_stats {

  executor = 'local'
  cpus = params.threads

  publishDir "${params.output_dir}/refinem", mode: "copy"

  input:
    path bins from Bins_filt_refinem_stats
    file ccs_filt_sort_bam
    file ccs_filt_sort_bam_index

  output:
    path "${params.sample_id}_scaffold_stats", type: "dir" into Refinem_scaffold_stats_out

  script:
    """
    ${REFINEM} \
    scaffold_stats \
    -x fa \
    --cpus ${params.threads} \
    --cov_all_reads \
    --cov_min_align ${params.min_align_len} \
    --cov_max_edit_dist ${params.cov_max_edit_dist} \
    ${params.contigs_fasta} \
    ${bins} \
    ${params.sample_id}_scaffold_stats \
    ${ccs_filt_sort_bam} \

    """
}

Refinem_scaffold_stats_out.into{Scafstats_for_outliers;
                                Scafstats_for_taxa}


// detect outlier scaffolds
// using refineM

process refinem_outliers {

  executor = 'local'
  cpus = 1

  publishDir "${params.output_dir}/refinem", mode: "copy"

  input:
    path stats from Scafstats_for_outliers

  output:
    path "${params.sample_id}_outliers", type: "dir" into Refinem_outliers_out

  script:
    """
    ${REFINEM} \
    outliers \
    --gc_perc ${params.gc_perc} \
    --td_perc ${params.td_perc} \
    --cov_perc ${params.cov_perc} \
    --individual_plots \
    --image_type svg \
    --point_size 36 \
    --label_font_size 12 \
    --tick_font_size 10 \
    --height 6 \
    --width 12 \
    ${stats}/scaffold_stats.tsv \
    ${params.sample_id}_outliers \

    """
}


// remove the outlier scaffolds
// with refinem filter_bins

process refinem_filter_outlier_bins {

  executor = 'local'
  cpus = 1

  publishDir "${params.output_dir}/refinem", mode: "copy"

  input:
    path bins from Bins_filt_filter_bins
    path outlier_bins from Refinem_outliers_out

  output:
    file "${params.sample_id}_filter_outlier_bins/bin_sizes.tsv" into bin_sizes
    path "${params.sample_id}_filter_outlier_bins", type: "dir" into Refinem_no_outliers

  script:
    """
    ${REFINEM} \
    filter_bins \
    -x fa \
    ${bins} \
    ${outlier_bins}/outliers.tsv \
    ${params.sample_id}_filter_outlier_bins &&
    { \
    echo -e "Bin\tSize(Mbp)"; \
    for i in `ls ${params.sample_id}_filter_outlier_bins/*fa`; \
    do \
    ${BIOAWK} \
    -c fastx \
    -v x=\${i} \
    '{sum += length(\$seq)} END {print x"\t"sum/1e6}' \${i}; \
    done | \
    sort -k2gr,2; } \
    > ${params.sample_id}_filter_outlier_bins/bin_sizes.tsv \

    """
}


// identify taxa contained in each MAG
// using the genome taxonomy database toolkit (gtdb-tk)

process run_gtdb_tk {

  executor = 'local'
  cpus = 1

  publishDir "${params.output_dir}/gtdb-tk", mode: "copy"

  input:
    path bins from Refinem_no_outliers

  output:
    path "${params.sample_id}_classify_wf" into Gtdb_out

  script:
    """
    ${GTDBTK} \
    classify_wf \
    --genome_dir ${bins} \
    --extension fa \
    --cpus ${params.threads} \
    --pplacer_cpus ${params.pplacer_threads} \
    --min_perc_aa ${params.min_perc_aa} \
    --min_af ${params.min_align_frac} \
    --out_dir ${params.sample_id}_classify_wf \
    --prefix ${params.sample_id} \

    """
}


// prepare taxonomy table

process taxonomy_table {

  executor = 'local'
  cpus = 1

  publishDir "${params.output_dir}/statistics", mode: "copy"

  input:
    file bin_sizes
    path Gtdb_out

  output:
    file "ALL.MAG_classification.tsv"

  script:
    """
    ${PYTHON3} \
    ${launchDir}/src/MAG-taxa-detection.py \
    --bac-summary ${Gtdb_out}/ALL.bac120.summary.tsv \
    --bin-sizes ${bin_sizes} \
    --output-file ALL.MAG_classification.tsv \

    """
}
