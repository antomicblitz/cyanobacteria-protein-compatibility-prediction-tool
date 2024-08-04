# cyanobacteria-compatibility-prediction-tool

## Overview

This tool analyzes and ranks proteins based on their potential suitability for expression in cyanobacteria. It evaluates proteins based on size, post-translational modifications (PTMs), cysteine content, and potential for disulfide bond formation.

## Requirements

- Python 3.6+
- Biopython library

## Installation

1. Clone this repository or download the script.
2. Ensure Python 3.6+ is installed on your system.
3. Install required dependencies:
pip install biopython
## Usage

1. Prepare an input CSV file named 'input_accessions.csv' with one protein accession number per line.
2. Run the script:
python protein_analysis.py

Two output files will be generated:
- 'ranked_proteins.csv': Contains proteins that passed screening, ranked by their composite score.
- 'excluded_proteins.csv': Lists proteins that did not pass screening and reasons for exclusion.

## Screening Process

The script screens proteins based on the following criteria:

1. Size: Proteins larger than 1500 amino acids are excluded.
2. Post-Translational Modifications (PTMs): Proteins with PTMs not typically possible in cyanobacteria are excluded or flagged.

### PTM Screening

- Excluded PTMs: glycosylation, ubiquitination, SUMOylation, complex acetylation, methylation, farnesylation, geranylgeranylation, palmitoylation, and myristoylation.
- Allowed PTMs: phosphorylation, disulfide bond formation.
- Proteins annotated with "Not [PTM]" or "Absence of [PTM]" are not excluded but may receive a lower confidence score.

## Scoring System

The composite score is calculated based on three factors:

1. Cysteine Density: (Number of cysteines / Protein length) * 100
2. Potential Disulfide Pairs: Count of cysteine pairs within 5-30 amino acids of each other
3. Spacing Score: Count of cysteine pairs spaced 5-30 amino acids apart

Each factor is normalized by dividing by the maximum value across all analyzed proteins. The composite score is then calculated as:
Composite Score = (0.3 * Normalized Cys Density) +
(0.4 * Normalized Potential Pairs) +
(0.3 * Normalized Spacing Score)
## Confidence Assignment

Confidence levels are assigned as follows:

- High: Passes all screening criteria without any concerns.
- Medium: Passes screening but has annotations suggesting absence of PTMs or other minor concerns.
- Low: Does not pass screening due to size or presence of excluded PTMs.

## Output Description

### ranked_proteins.csv

- Accession: Protein accession number
- Name: Protein name/description
- Length: Number of amino acids
- Cys Density: Cysteines per 100 amino acids
- Potential Pairs: Number of potential disulfide bond-forming cysteine pairs
- Spacing Score: Score based on cysteine spacing
- Composite Score: Final ranking score
- Screening Result: "Passed" for all entries in this file
- Confidence: High, Medium, or Low

### excluded_proteins.csv

- Accession: Protein accession number
- Name: Protein name/description
- Exclusion Reason: Specific reason for exclusion
- Confidence: Always "Low" for excluded proteins

## Limitations and Considerations

- This tool provides theoretical predictions and should be used as a starting point for selecting proteins for experimental validation.
- The scoring system prioritizes cysteine content and potential disulfide bonds, which may not be the only factors determining successful expression in cyanobacteria.
- PTM annotations in protein databases may not be exhaustive or may contain errors.
- While the tool considers some cyanobacteria-specific PTM capabilities, the exact PTM landscape can vary between different cyanobacterial species.
- Experimental validation is crucial to confirm actual expression success and proper protein folding in cyanobacteria.
- Other factors such as codon usage, protein solubility, and specific folding requirements are not considered in this analysis.

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
