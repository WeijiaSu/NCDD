# NCDD - Nanopore Circle DNA Detector

The NCDD.py script is designed to analyze nanopore sequencing data to identify and classify circular DNA. The script analyzes alignment files and extracts information about circular DNA based on user-defined criteria.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Usage](#usage)
- [Function Descriptions](#function-descriptions)
- [Output Files](#output-files)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites

- **Python 3.x**
- **Libraries and packages**:
  - os
  - argparse
  - re
  - Bio (from BioPython)
  - pysam
  - numpy
  - warnings
  - cigar (optional)

## Usage

```bash
python NCDD.py -bam [BAMFILE] -js [JUNSIZE] -Prefix [PREFIX] -ref [REFERENCE]

## Parameters:

-bam or --bamFile: Input BAM file. (required)
-js or --JunSize: Size of the junction. (default is 100)
-Prefix or --Prefix: Prefix for output files. (default is None)
-ref or --reference: Reference genome in FASTA format. (required)

## Function Descriptions
map_ratio: Calculates the ratio of mapped regions in the read.
getChimeric_reads: Extracts chimeric reads based on mapping ratios.
FilterReads: Filters out linearly aligned reads and retains potential circular DNA reads.
Junction: Checks if two aligned fragments form a valid junction.
isCircle: Determines if the sequence forms a circular DNA based on alignment coordinates.
CompleteCopy: Determines if a sequence has a complete copy of the potential circular region.
CircleType: Classifies the type of circle based on alignment patterns.
getCircle: Identifies circular reads from the data.
GetReads: Main function that orchestrates the entire analysis.


## Output Files
PREFIX_rc90.tsv: File containing reads that have more than 90% alignment coverage.
PREFIX_LiAg.tsv: File containing linear alignments.
PREFIX.candi.tsv: File containing candidate reads for circular DNA.
PREFIX_circles.txt: File containing reads classified as circles.
PREFIX_candi.tsv_circleAnalyze.txt: Analysis of potential circular reads.
Contributing
For any suggestions, bug reports, or any other feedback, please open an issue on GitHub.

