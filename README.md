# NCDD.py - Nanopore Circle DNA Detector

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

