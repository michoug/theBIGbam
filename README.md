# Installation

## Basic installation

## Additional dependencies regarding mapping
If you want to perform mapping using the internal scripts provided, you will need several tools installed on your path:
- minimap2 and samtools for read mapping commands
- seqkit and either pharokka or bakta for assembly annotation command
An easy way is to use a conda environment

mamba env create -f mgfeatureviewer_env.yaml
conda activate mgfeatureviewer

git clone https://github.com/bhagavadgitadu22/MGFeatureViewer
cd MGFeatureViewer
python -m pip install -e .

mgfeatureviewer -h
# or alternatively: python -m mgfeatureviewer.cli -h

rust install?

# Usage
If you do not have an assembly annotation file (.genbank) or sorted mapping files with MD tags (.bam), you can start with the scripts provided in the "Assembly annotation" and "Mapping" sections below to obtain the required input files.

## Database computation
This is the core of the MGFeatureViewer pipeline. It takes a genbank file containing your contigs of interest and a list of bam files containing the mappings against all or some of your contigs of interest. You should specify the --circular flag if you did a circular mapping as explained in the "Mapping" section. Only contigs where more than --min_coverage % of their length are covered by reads are considered during calculations to speed up the program. Moreover, you can balance between sensitivity and compactness with the --compress_ratio option: the higher the lighter the final database but with a risk of missing positions with potentially relevant features.

Fast Rust calculations are performed on your bam files to extract the relevant values depending on the features you asked for. 3 kinds of features can be requested at the moment:
- coverage: computes the coverage of primary, secondary and supplementary reads per position
- phagetermini: looks at the starts and ends of mapped reads, this is useful in particular to identify the termini of phages (see the [PhageTerm publication](https://www.nature.com/articles/s41598-017-07910-5)). A stringent coverage is computed with only reads starting with a match (and ending with a match in the case of long reads). The number of reads starting and ending at each position is computed as well as the ratio tau = (reads_starts + reads_ends) / stringent_coverage
- assemblycheck: computes a number of useful metrics to assess the completion and contamination of a contig. Clippings, indels and mismatches are calculated. Read lengths are computed for long reads, while for paired short reads you compute insert sizes and wrong orientations between a pair of reads

For each sample, Rust calculations are computed for all contigs considered present in the bam file (more than --min_coverage % of their length are covered by reads). Further reduction of the final database is performed post-calculation during a compressing step. For the different types of coverage, compressing is done via a Run-Length Encoding (RLE logic): consecutive positions with similar values (within the --compress_ratio threshold) are grouped together as one single entry. This allows to keep the overall shape of the signal while drastically reducing the number of entries to store. For other features such as phage termini and assembly check metrics, values are compared to the local coverage: if they are below the --compress_ratio threshold relative to the coverage at that position, they are considered insignificant and filtered out. This allows to keep only relevant peaks in the data while discarding noise. If consecutive positions are conserved, RLE logic is applied to group them together. 

The output is an SQLite database containing all those values. This database is typically 1,000 times smaller than the original bam files while maintaining the essential characteristics of the mapping data.

Example of command to compute the database:
```sh
# rm examples/outputs/HK97/HK97.db
# maturin develop

DIR="examples/inputs/HK97"
mgfeatureviewer calculate -t 4 -g ${DIR}/HK97_GCF_000848825.1_pharokka.gbk -a pharokka -b ${DIR}/bams -m coverage,phagetermini,assemblycheck -c -o examples/outputs/HK97/HK97.db --circular --compress-ratio=10

# alternatively you can compute the database directly in Rust with
cargo run -- -t 4 -g ${DIR}/HK97_GCF_000848825.1_pharokka.gbk -a pharokka -b ${DIR}/bams -m coverage,phagetermini,assemblycheck -o examples/outputs/HK97/HK97.db --circular --compress-ratio=10
```

## Visualisation
Example command:
```sh
cp examples/outputs/HK97/HK97.db ~/HK97.db
mgfeatureviewer serve --db ~/HK97.db --port 5006
```

## Updating the database
Example command:
```sh
mgfeatureviewer add-variable examples/outputs/HK97/HK97.db test bars "#f60820" "Test Bar" examples/inputs/HK97/variable_test.csv
mgfeatureviewer add-variable examples/outputs/HK97/HK97.db test2 curve "#aef1c2" "Test Curve" examples/inputs/HK97/variable_test2.csv
```

# Utils

## Mapping
You can start the pipeline with your reference contigs and read files. You can map one sample only or several samples at once. Sorted ".bam" files are generated for each file along an index file ".bam.bai".

To map several samples, you need to provide a csv file containing 4 columns separated by commas: "read1,read2,sequencing_type,assembly_file". Only read1 is compulsory: read2 is only used for paired-reads, while sequencing_type and assembly_file can be provided via the global options -s and -a whether their values should be the same for all samples.

Example of csv file (available in examples/inputs/HK97/mapping_rows.csv):
```
examples/inputs/HK97/HK97_R1_illumina.fastq.gz,examples/inputs/HK97/HK97_R2_illumina.fastq.gz,short,examples/inputs/HK97/HK97_GCF_000848825.1.fasta
examples/inputs/HK97/HK97_nanopore.fastq.gz,,long,examples/inputs/HK97/HK97_GCF_000848825.1.fasta
```

The mapping proposed here allows a circularized mapping via the --circular flag. All contigs are concatenated to themselves, doubling the size of the assembly. It allows reads to map at the ends of the contigs without breaking. During the rest of the pipeline all calculation take into account that the mapping was circular if the --circular flag is provided: mappings at positions beyond the end of a contig are counted at that position minus the size of the contig instead.

Another particularity of this mapping is that MD tags are computed using samtools. MD tags are useful to quickly identify mismatches among the mapped reads.

Examples of commands to map a single sample or several:
```sh
DIR="examples/inputs/HK97"
mgfeatureviewer mapping-per-sample -s "short" -r1 "${DIR}/HK97_R1_illumina.fastq.gz" -r2 "${DIR}/HK97_R2_illumina.fastq.gz" -a "${DIR}/HK97_GCF_000848825.1.fasta" -o "${DIR}/HK97_GCF_000848825.1_with_MD.bam" --circular

mgfeatureviewer mapping-all-samples --csv "${DIR}/mapping_rows.csv" -o "${DIR}/bams" --circular
```

## Assembly annotation
TODO 




mgfeatureviewer plot-per-sample -d examples/outputs/HK97/HK97.db -v "Coverage,Phage termini,Assembly check,test" --contig NC_002167.1 --sample HK97_GCF_000848825.1_with_MD --html examples/outputs/HK97/HK97_illumina_per_sample.html
mgfeatureviewer plot-all-samples -d examples/outputs/HK97/HK97.db -v "Coverage" --contig NC_002167.1 --html examples/outputs/HK97/HK97_illumina_all_samples.html



# How to interpret the plots?

## On the coverage:
coverage and coverage_reduced use different kind of filtering.

coverage:
- only counts one cover if a read truly matches the given basepair
coverage_reduced:
- considers less reads than coverage because it only considers reads starting and finishing with a match (not clipping authorised at the ends of the read)
- counts one cover for all the positions between the first basepair and the last basepair of the read

Thus for a given position both coverage can be higher than the other one depending on the scenarii:
- coverage > coverage_reduced when part of the genome is missing in the assembly resulting in a lot of clippings at the position considered
- coverage_reduced > coverage when the basepair is missing in some of the organism' population

A particularly low coverage_reduced compared to coverage suggests your reads were not trimmed properly before mappings: you likely still have adapter sequences.


# TODO:
explain what are different features
explain database structure
explain how to interrogate database with sql
assembly annotation explanation

mention how secondary is an approximation in case of doubled mapping because some primary at the ends of the contig would not have associated secondary. as i keep a minimum of 0 secondary, i will not see if some true secondary appear there. but this is minor and allows to go way faster than checking each read individually to remove secondary associated exactly to a primary read

tab to complete in plot  Contig and Sample if only one value possible


need to rethink avout what i plot
secondary and supplemnentary if i had them for minimap2
i get insert_sizes and bad_orientations messed up in circular mode
also think about proper way to calculate sequencing_type
and parameters i should pass to assembly_annotation