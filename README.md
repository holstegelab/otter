# *otter*
A generic targeted local-assembler for regions of interest for (aligned) long-read data. It automaticaly adapts to the local sequencing error-rates and coverages, and hence can be used on data from both PacBio and ONT instruments. It employs a genotyper to enable comparison of the assembled sequences across multiple samples. 

It was initially devoloped for genotyping tandem repeats across one or more genomes, but can also be used for other structural variants such as transposable elements.

If given whole-genome alignments, *otter* can be used to compare and genotype regions of interest across multiple *de novo* assemblies.

## Dependencies

Required:
* A c++17 supporting compiler

Included during install:
* [WFA2-lib](https://github.com/smarco/WFA2-lib)
* [hclust-cpp](https://github.com/cdalitz/hclust-cpp/tree/master)

Optional, required for downstream analysis (VCF-generation and joint-genotyping):
* [samtools](https://www.htslib.org/download/)

## Installation

```bash
git clone --recursive https://github.com/holstegelab/otter.git && cd otter
# build WFA2-lib dependency
cd include/WFA2-lib; make clean setup lib_wfa; cd ../../
mkdir build
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
otter assemble -b regions.bed -R sample_name --fasta reads.bam > assembly.fa
```

If an indexed reference genome is provided (*ref.fa*), *otter* will perform local read re-alignment on-the-fly if it detects potential suboptimal alignments due to highly divergent sequence between the region of interest in the read and a reference:
```bash
otter assemble -b regions.bed -r ref.fa -R sample_name --fasta reads.bam > assembly.fa
```

You can view the corresponding sequences of individual reads per region by using the ```--reads-only``` parameter.

Although the above examples output in FASTA-format, it is recommended to output in SAM/BAM to enable VCF generation and joint-genotyping (see section below):

```bash
#note: samtools is required
otter assemble -b regions.bed -r ref.fa -R sample_name reads.bam | samtools view -bh | samtools sort > assembly.bam
samtools index assembly.bam
#generate VCF file
otter genotype -b regions.bed -r ref.fa assembly.bam > assembly.vcf
```

### Joint-genotyping

*otter* can perform joint-genotyping of target-regions across multiple samples (e.g. a cohort) to enable uniform representations of structural variants. You can do this in three steps:

1. Locally assemble regions of interest with *otter* in BAM-format (see above).

2. Generate a single merged-BAM:
```bash
samtools merge -pco cohort.bam sample1.otter.bam sample2.otter.bam sample3.otter.bam ...
samtools index cohort.bam
```
3. Perform joint-genotyping and output in VCF-format:
```bash
otter genotype -b regions.bed -r ref.fa merged.bam > merged.vcf
```

### Genotyping *de novo* assemblies

It is possible to use the same workflow above to genotype regions of interest in one or more aligned whole-genome *de novo* assemblies. The whole-genome alignments must be in BAM-format. To genotype:

```bash
#for FASTA-format
otter wgat -b regions.bed --fasta aligned_assembly.bam > aligned_assembly.fa
#for BAM-format
otter wgat -b regions.bed aligned_assembly.bam | samtools view -bh | samtools sort > aligned_assembly.bam
```
### Metadata
There is metadata available for each allele-sequence in the auxillarary tags (BAM or FASTA):
| Tag | Type | Description |
| --- | ---- | ----------- |
| RG | Z | Sample name |
| ta | Z | Region of interest (BED-format) |
| ic | i | Initial set of unique allele-sequences during assembly |
| tc | i | Total reads in region of interest |
| sc | i | Total spanning reads in this region |
| ac | i | Total reads assigned to this unique allele-sequence |
| se | f | Standard-error of pairwise sequence dissimilarity for all reads contributing to this unique allele-sequence |

Additional tags are available for alignment-based output (e.g. ```wgat``` command or ```--reads-only``` parameter in ```assemble``` command):
| Tag | Type | Description |
| --- | ---- | ----------- |
| sp | A | Alignment spanning status |

Definitions:
| Value | Description |
| --- | ----------- |
| b | Alignment spans from both flanks |
| l | Alignment spans only from left-flank |
| r | Alignment spans only from right-flank |
| n | Alignment does not span from either flank |

### Limitations
* *otter* was developed for characterising structural variation, and hence it may not be able distinguish highly similar sequences (e.g. few nucleotide differences). In these cases, it will collapse the assembly into a single allele-sequence.
* The backbone of the assembly relies on spanning-reads, and hence there is a size limitation (spannable by the long-reads).

### Citations
[See our preprint (update coming soon)](https://www.biorxiv.org/content/10.1101/2024.03.15.585288v1)