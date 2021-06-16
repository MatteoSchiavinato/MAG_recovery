# Workflow for the annotation of metagenomic assembled genomes (MAGs) from pacbio read data

This workflow takes as input a set of pacbio hifi reads and a set of contigs assembled from them, and produces as output a taxonomically annotated set of contig clusters. There are many steps required to filter and cluster the contigs properly, which are performed with the [metabat2](https://doi.org/10.7717/peerj.7359), [checkm](https://ecogenomics.github.io/CheckM/) and [refinem](https://github.com/dparks1134/RefineM) programs. The taxonomic annotation is performed with [gtdb-tk](https://github.com/Ecogenomics/GTDBTk). The clusters represent the closest representation of individual genomes in the sample. These clusters are called metagenomic assembled genomes (MAGs).

### Input requirements

This workflow only takes two input sources:
- The assembled contigs in FASTA format
- The raw reads in FASTA format

### Dependencies

The pipeline depends on the following tools:

| Software    | Version | Type        | Link                                    |
|-------------|---------|-------------|-----------------------------------------|
| Nextflow    | 21.04.1 | Interpreter | https://www.nextflow.io/                |
| Minimap2    | 2.18    | Program     | https://github.com/lh3/minimap2         |
| Samtools    | 1.12    | Program     | http://www.htslib.org/download/         |
| Metabat2    | 2.15    | Program     | https://anaconda.org/bioconda/metabat2  |
| Checkm      | 1.1.3   | Program     | https://ecogenomics.github.io/CheckM/   |
| Refinem     | 0.1.2   | Program     | https://github.com/dparks1134/RefineM   |
| Gtdb-tk     | 1.5.0   | Program     | https://github.com/Ecogenomics/GTDBTk   |
| Python      | 3.8     | Interpreter | https://www.python.org/downloads/       |
| Pandas      | 1.2.4   | Module      | https://pandas.pydata.org/              |


### Run the pipeline

To see the pipeline help section, run:

```
nextflow run main.nf --help
```

Then, run the pipeline with a command like this one:

```
nextflow \
run \
main.nf \
--output_dir /path/to/output \
--threads 48 \
--contigs_fasta /path/to/contigs.fasta \
--ccs_reads /path/to/ccs_reads.fasta \
--sample_id <sample_id> \

```

##### Nextflow tips

- You can change all parameters from command line or from the `nextflow.config` file. Not all parameters are shown with `--help`, so check the `nextflow.config` file too.
- You can resume a crashed run by re-running the same command and adding `-resume` as an option after `run`. More info [here](https://www.nextflow.io/docs/latest/getstarted.html).
- You can specify where to save the temporary files of the pipeline by specifying a `-work-dir` directory right after `run`. More info on that and on other options available in `nextflow run` can be found [here](https://www.nextflow.io/docs/latest/cli.html#clean).
- The `work` directory tends to become quite crowded and full thousands of internal nextflow files, so every now and then make sure you clean it with `nextflow clean`. More info [here](https://www.nextflow.io/docs/latest/cli.html#clean)

### Steps and output

The pipeline is articulated in several steps, which are described below.

##### Contig binning

In the first step, the contigs passed via command line with the `--contigs_fasta` argument are passed to **metabat2**. A series of parameters can be set to control this step (see `--help`). The most important is perhaps the minimum cluster total length (`--min_cluster_tot_len`). By default, this parameter is set to "100000" (0.1 Mbp). Considering what MAGs are, you perhaps don't want to retain clusters that are too small. You may want to consider the default value as the lowest possible value, and if you feel like it, increase it.

The output of this step is contained in the `metabat2` folder inside the `--output_dir`. The folder will contain a subfolder named like the sample passed with `--sample_id`, plus the `_out` suffix. The sample ID folder will contain all the bins in FASTA format, with an associated number.

##### Bin checking

The bins are then passed by the pipeline to the **checkm** program, which performs a series of analysis on the bins to determine their completeness, heterogeneity, and contamination. It also attempts to place them in the tree of life by looking at a series of marker genes. However, since the pipeline will then perform a more thorough bin tree placing, the taxonomic results of this step are not used by this pipeline. Only the contamination, heterogeneity and completeness are used to filter the bins.

*Completeness*: based on how many marker genes are found from the `checkm lineage_wf` program.
*Contamination*: based on how many of these marker genes are found *twice* or more.
*Heterogeneity*: based on sequence identity between found marker genes.

See the [checkm manual](https://ecogenomics.github.io/CheckM/) for details, this wasn't something created by me so it's better to look at the manual of the original program.

The output of this step is contained in the `checkm` folder inside the `--output_dir`. The folder will contain two subfolders named like the sample passed with `--sample_id`, plus the `_raw_lineage` and the `_filtered_lineage` suffixes. Both folders contain intermediate files of the `checkm lineage_wf` program.

The `checkm` directory will also contain the `checkm_lineage_wf_results.tsv` which is the real result of a filtering step performed with a custom python script that you can find in the `/src` subdirectory of this workflow. This file contains the **completeness**, **contamination** and **heterogeneity** values for each bin that was retained after filtering. The filtering parameters can be passed with the `--max_contamination`, `--min_completeness`, and `--max_heterogeneity` options. As mentioned before, the taxonomic annotation that you'll find in this file for each bin (2nd column) won't be used because a more thorough taxonomic annotation will be performed later on in the pipeline.

##### Bin refining

At this stage, the pipeline refines the contigs contained in each bin based on their [tetranucleotide frequencies](https://doi.org/10.1002/elps.1150190412), coverage from ccs reads, and GC content distribution. First, the pipeline maps the original ccs reads with **minimap2** against the contigs. Then it extracts a filtered and sorted bam file with `samtools view -b -h -F 0x0100 -F 0x4`, removing secondary alignments and unmapped records. This file is then passed to `refinem scaffold_stats` together with the bins produced by `checkm lineage_wf`. The "scaffold stats" tool analyzes the frequency of every tetranucleotide contained in the contigs of each bin. It also uses the bam file to produce coverage profiles per position. It also calculates the GC content distribution. These three informations are then used to detect whether any contig in any bin seems to be an outlier compared to the other contigs in its bin. This operation is performed with `refinem outliers` (see [Parks et al., 2017](https://doi.org/10.1038/s41564-017-0012-7) for details). The outlier bins are removed with `refinem filter_bins`.

The tetranucleotide frequency percentile, the GC content percentile, and the coverage percentile can be controlled with command-line parameters: `--td_perc`, `--gc_perc`, `--cov_perc`. The default values are 95, 95, 50, respectively. The higher the value (1-100) the more the chances that a contig is considered an outlier. The original paper from this method suggests to use a 95-th percentile threshold for GC content and tetranucleotide frequency, as these two properties of a genome are somewhat stable and species-specific.

The output of this step is contained in the `refinem` folder inside the `--output_dir`. The folder will contain three subfolders named like the sample passed with `--sample_id`, plus the `_scaffold_stats`, the `_outliers` and the `_filter_outlier_bins` suffixes.

The important output of this step is contained in the `*_filter_outlier_bins` folder. This folder will contain each bin that survived the filtering in FASTA format. These bins are the input of the next step.

##### Bin taxonomic annotation

In the end, a taxonomic annotation of the bins is performed with [gtdb-tk](https://github.com/Ecogenomics/GTDBTk). This tool makes use of the most modern resources for metagenomic taxonomic annotation. For further information, you can read the following papers: [Parks et al., 2018](http://dx.doi.org/10.1038/nbt.4229), [Parks et al., 2020](https://doi.org/10.1038/s41587-020-0501-8).

The results are contained in the `gtdb-tk` and in the `statistics` folders inside the `--output_dir`. The `gtdb-tk` folder contains all the files produced by **gtdb-tk**. However, the real result is in `statistics` and is the `ALL.MAG_classification.tsv` file. This file contains all the known taxonomic annotation for each MAG, including its size.
