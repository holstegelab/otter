# *otter*
A generic targeted local-assembler for regions of interest for (aligned) long-read data. It automaticaly adapts to the local sequencing error-rates and coverages, and hence can be used on data from both PacBio and ONT instruments. The assembled sequences from one or more samples can be automatically stored in SAM/BAM-format, enabling compatibility with downstream SAM/BAM-compatible tools. It employs a zygosity-aware genotyper to enable comparison of the assembled sequences across multiple samples. 

It was initially devoloped for genotyping tandem repeats across one or more genomes, but can also be used for other structural variants such as transposable elements.

If given whole-genome alignments, *otter* can be used to compare and genotype regions of interest across multiple *de novo* assemblies.

## Dependencies

* A c++17 supporting compiler
* [htslib-1.3+](https://www.htslib.org/download/)
* [WFA2-lib](https://github.com/smarco/WFA2-lib) (included during install)
* [hclust-cpp](https://github.com/cdalitz/hclust-cpp/tree/master) (included)

Note that *htslib* needs be available in your library paths. You can do this by adding the *htslib* path in your *.bashrc* or *.bash_profile* file:

```bash
export CPATH=$CPATH:<your_local_htslib_path>
export LIBRARY_PATH=$LIBRARY_PATH:<your_local_htslib_path>
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<your_local_htslib_path>
```

## Installation

```bash
git clone https://github.com/holstegelab/otter.git && cd otter
mkdir build
make packages
make
```

## Usage

See commands and version by running:
```bash
./build/otter
```

### Targeted local assembly

Given a BED-formatted file (*regions.bed*) and aligned long-reads (*reads.bam*), perform targeted local assembly for each region of interest, and will output all unique allele-sequences:
```bash
otter assemble -b regions.bed reads.bam > assembly.fa
```
If an indexed reference genome is provided (*ref.fa*), *otter* will perform local read re-alignment on the fly if it detects potential suboptimal alignments due to highly divergent sequence between the region of interest in the read and a reference:
```bash
otter assemble -b regions.bed -r ref.fa reads.bam > assembly.fa
```

You can view the corresponding sequences of individual reads per region by using the ```--reads-only``` parameter.

By default, the output is FASTA-formatted with the following format in each sequence name: \<BED-string of region\> \<allele-coverage\> \<total coverage\> \<standard error\>

We recommend outputting in SAM-format for enabling joint multi-sample genotyping and comparison (see section below). A unique sample-name is required:
```bash
otter assemble -b regions.bed -r ref.fa --sam -R sample_name > assembly.sam
```
Note that you can automatically compress and sort to BAM-format with [*samtools*](https://github.com/samtools/samtools), enabling low-storage requirement and downstream usability with other SAM/BAM-compatible tools:
```bash
otter assemble -b regions.bed -r ref.fa --sam -R sample_name reads.bam | samtools view -bh | samtools sort > assembly.bam
samtools index assembly.bam
```
There is additional meta-data in each unique allele-sequence in the auxillarary tags:
| Tag | Type | Description |
| --- | ---- | ----------- |
| RG | Z | Sample name |
| ta | Z | Region of interest (BED-format) |
| ic | i | Initial set of unique allele-sequences during assembly |
| ac | i | Total reads assigned to this unique allele-sequence |
| tc | i | Total reads in region of interest |
| se | f | Standard-error of pairwise sequence dissimilarity for all reads contributing to this unique allele-sequence |

### Genotyping *de novo* assemblies

It is possible to use the same workflow above to genotype regions of interest in one or more (whole-genome) *de novo* assemblies. Given whole-genome alignments in BAM-format, *otter* will extract the corresponding allele-sequences by using the ```--wga``` parameter.

### Joint multi-sample genotyping

You can merge separate BAM-formatted *otter* assemblies of multiple samples into a single BAM-file:
```bash
samtools merge -o merged.bam sample1.bam sample2.bam sample3.bam ...
samtools index merged.bam
```

With the merged BAM-file, you can perform multi-sample genotyping and comparison of the local assemblies per region. We implement a few strategies in *otter*:

#### Genotyping by clustering 

For each region, *otter* will perform zygosity-aware pairwise-sequence alignment, calculate sequence dissimilarity, and perform hierchical clustering up-to a given maximum dissimilarity threshold, *e<sub>max</sub>* (by default 0.05). Each cluster is given an index (e.g. 0, 1), and used to assign the genotype for a sample (e.g. 0/0, 0/1, 1/1):
```bash
otter genotype -b regions.bed -e e_max merged.bam > genotypes.tsv
```
The output is a TSV-formatted file: \<BED-string of region\> \<sample_name\> \<Genotype\>

Note that highly variable genomic regions (e.g. tandem repeats) may be multi-allelic leading to >2 unique allele-sequences.

#### Genotyping by length

You can output zygosity-aware sequence-lengths per-sample per-region with the same output format as above but replacing the \<Genotype\> column with three separate columns: \<min length\> \<max length\> \<sum length\>
```bash
otter genotype -b regions.bed -l merged.bam > lengths.tsv
```

### Limitations
* *otter* was developed for characterising structural variation, and hence it may not be able distinguish highly similar sequences (e.g. few nucleotide differences). In these cases, it will collapse the assembly into a single allele-sequence.
* The backbone of the assembly relies on spanning-reads, and hence there is a size limitation (spannable by the long-reads).

### Citations
Pre-print coming soon.