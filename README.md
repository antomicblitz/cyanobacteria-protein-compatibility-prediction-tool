# Cyanobacteria Protein Compatibility Prediction Tool

## Overview

This tool analyzes and ranks proteins based on their potential suitability for expression in cyanobacteria. It evaluates proteins based on size, specific post-translational modifications (PTMs), and cysteine content. The tool provides detailed information about why each protein passes or fails the screening process.

## Requirements

- Python 3.6+
- Biopython library

## Installation

1. Clone this repository or download the script.
2. Ensure Python 3.6+ is installed on your system.
3. The script will create a virtual environment and install the required dependencies automatically.

## Usage

1. Prepare an input CSV file named 'input_accessions.csv' with one protein accession number per line.
2. Run the script:
python protein_analysis.py
Copy
Two output files will be generated:
- 'ranked_proteins.csv': Contains proteins that passed screening, ranked by their cysteine density.
- 'excluded_proteins.csv': Lists proteins that did not pass screening and reasons for exclusion.

## Screening Process

The script screens proteins based on the following criteria:

1. Size: Proteins larger than 1500 amino acids are excluded.
2. Post-Translational Modifications (PTMs): Proteins with PTMs not typically possible in cyanobacteria are excluded.

### PTM Screening

- Excluded PTMs: complex glycosylation, ubiquitination, SUMOylation, farnesylation, geranylgeranylation, palmitoylation, and myristoylation.
- Compatible PTMs: phosphorylation, simple glycosylation, methylation, acetylation, and disulfide bonds.
- Proteins annotated with "Not [PTM]" or "Absence of [PTM]" for excluded PTMs are not excluded.

## Ranking System

Proteins are ranked based on their cysteine density, calculated as:
Cysteine Density = (Number of cysteines / Protein length) * 100

## Confidence Assignment

Confidence levels are assigned as follows:

- High: Passes all screening criteria without any concerns.
- Low: Does not pass screening due to size or presence of excluded PTMs.

## Output Description

### ranked_proteins.csv

- Accession: Protein accession number
- Name: Protein name/description
- Length: Number of amino acids
- Cys Density: Cysteines per 100 amino acids
- Screening Result: "Passed" for all entries in this file
- Confidence: High or Low
- Pass Reason: Detailed explanation of why the protein passed, including compatible PTMs found and high cysteine density if applicable

### excluded_proteins.csv

- Accession: Protein accession number
- Name: Protein name/description
- Exclusion Reason: Specific reason for exclusion (size or incompatible PTM)
- Confidence: Always "Low" for excluded proteins

## Interpreting Results

- Proteins with high cysteine density (>5%) are noted in the Pass Reason column.
- Compatible PTMs found in the protein are listed in the Pass Reason column.
- Proteins without detected compatible PTMs or high cysteine density will be marked as "No incompatible PTMs found".

## Limitations and Considerations

- This tool provides theoretical predictions and should be used as a starting point for selecting proteins for experimental validation.
- The detection of PTMs is based on annotations in the protein database and may not be exhaustive.
- While the tool considers some cyanobacteria-specific PTM capabilities, the exact PTM landscape can vary between different cyanobacterial species.
- Experimental validation is crucial to confirm actual expression success and proper protein folding in cyanobacteria.
- Other factors such as codon usage, protein solubility, and specific folding requirements are not considered in this analysis.

## Customization

You can modify the following parameters in the script:

- `max_size`: Maximum allowed protein size (default: 1500 amino acids)
- `excluded_ptms`: List of PTMs incompatible with cyanobacteria
- `compatible_ptms`: Dictionary of PTMs compatible with cyanobacteria
- Cysteine density threshold for "high cysteine density" classification (default: 5%)

## Troubleshooting

If you encounter any issues:

1. Ensure you have internet access, as the script fetches protein information from NCBI.
2. Check that your input CSV file is correctly formatted with one accession number per line.
3. Verify that you have the necessary permissions to create files in the script's directory.

## License

This project is licensed under the MIT License:

MIT License

Copyright (c) 2024 Antonio Lamb

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Contact

antoniolamb at gmail dot com
