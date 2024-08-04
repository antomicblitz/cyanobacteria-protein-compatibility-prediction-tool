import csv
import os
import venv
import subprocess
import sys
from Bio import Entrez, SeqIO
import re

def create_virtual_environment():
    venv_dir = 'protein_analysis_env'
    venv.create(venv_dir, with_pip=True)
    
    if sys.platform == 'win32':
        python_path = os.path.join(venv_dir, 'Scripts', 'python.exe')
        pip_path = os.path.join(venv_dir, 'Scripts', 'pip.exe')
    else:
        python_path = os.path.join(venv_dir, 'bin', 'python')
        pip_path = os.path.join(venv_dir, 'bin', 'pip')
    
    subprocess.run([pip_path, 'install', 'biopython'])
    
    print(f"Virtual environment created at {venv_dir}")
    return python_path

def run_in_virtual_environment(python_path):
    script_path = os.path.abspath(__file__)
    subprocess.run([python_path, script_path, '--run-analysis'])

def fetch_protein_info(accession):
    Entrez.email = "your.email@example.com"  # Replace with your email
    handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    return record.description, str(record.seq), record.features

def analyze_protein(sequence, features, max_size=1500):
    if len(sequence) > max_size:
        return None, f"Protein too large ({len(sequence)} aa)", "Low"

    excluded_ptms = [
        "glycosyl", "ubiquitin", "sumo", "acetyl", "methyl",
        "farnesyl", "geranylgeranyl", "palmitoyl", "myristoyl"
    ]
    
    allowed_ptms = ["phospho", "disulfide"]
    
    ptms = [feature for feature in features if feature.type == "Site"]
    confidence = "High"
    for ptm in ptms:
        ptm_description = ptm.qualifiers.get('note', [''])[0].lower()
        if any(allowed_ptm in ptm_description for allowed_ptm in allowed_ptms):
            continue
        if any(re.search(rf'\b{excluded_ptm}\w*', ptm_description) for excluded_ptm in excluded_ptms):
            if "not " in ptm_description or "absence of " in ptm_description:
                confidence = "Medium"
                continue
            return None, f"Excluded PTM found: {ptm_description}", "Low"
    
    cys_count = sequence.count('C')
    cys_density = (cys_count / len(sequence)) * 100
    
    cys_positions = [i for i, aa in enumerate(sequence) if aa == 'C']
    potential_pairs = sum(1 for i, pos in enumerate(cys_positions[:-1])
                          if 5 <= cys_positions[i+1] - pos <= 30)
    
    spacing = {}
    for i in range(len(cys_positions) - 1):
        space = cys_positions[i+1] - cys_positions[i]
        spacing[space] = spacing.get(space, 0) + 1
    
    spacing_score = sum(count for space, count in spacing.items() if 5 <= space <= 30)
    
    return {
        'length': len(sequence),
        'cys_density': cys_density,
        'potential_pairs': potential_pairs,
        'spacing_score': spacing_score
    }, "Passed", confidence

def rank_proteins(accessions):
    results = []
    excluded = []
    for accession in accessions:
        try:
            name, sequence, features = fetch_protein_info(accession)
            analysis, message, confidence = analyze_protein(sequence, features)
            if analysis:
                analysis['accession'] = accession
                analysis['name'] = name
                analysis['screening_result'] = message
                analysis['confidence'] = confidence
                results.append(analysis)
            else:
                excluded.append({
                    'accession': accession,
                    'name': name,
                    'exclusion_reason': message,
                    'confidence': confidence
                })
                print(f"Excluded {accession}: {message} (Confidence: {confidence})")
        except Exception as e:
            excluded.append({
                'accession': accession,
                'name': 'Error',
                'exclusion_reason': str(e),
                'confidence': 'Low'
            })
            print(f"Error processing {accession}: {str(e)}")
    
    if not results:
        print("No proteins passed the screening criteria.")
        return [], excluded

    max_density = max(r['cys_density'] for r in results)
    max_pairs = max(r['potential_pairs'] for r in results)
    max_spacing = max(r['spacing_score'] for r in results)
    
    weights = {
        'cys_density': 0.3,
        'potential_pairs': 0.4,
        'spacing_score': 0.3
    }
    
    for r in results:
        r['composite_score'] = (
            weights['cys_density'] * (r['cys_density'] / max_density) +
            weights['potential_pairs'] * (r['potential_pairs'] / max_pairs) +
            weights['spacing_score'] * (r['spacing_score'] / max_spacing)
        )
    
    return sorted(results, key=lambda x: x['composite_score'], reverse=True), excluded

def main():
    input_file = 'input_accessions.csv'
    output_file = 'ranked_proteins.csv'
    excluded_file = 'excluded_proteins.csv'

    with open(input_file, 'r') as f:
        reader = csv.reader(f)
        accessions = [row[0] for row in reader]

    ranked_proteins, excluded_proteins = rank_proteins(accessions)

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Accession', 'Name', 'Length', 'Cys Density', 'Potential Pairs', 'Spacing Score', 'Composite Score', 'Screening Result', 'Confidence'])
        for protein in ranked_proteins:
            writer.writerow([
                protein['accession'],
                protein['name'],
                protein['length'],
                f"{protein['cys_density']:.2f}",
                protein['potential_pairs'],
                protein['spacing_score'],
                f"{protein['composite_score']:.4f}",
                protein['screening_result'],
                protein['confidence']
            ])

    with open(excluded_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Accession', 'Name', 'Exclusion Reason', 'Confidence'])
        for protein in excluded_proteins:
            writer.writerow([
                protein['accession'],
                protein['name'],
                protein['exclusion_reason'],
                protein['confidence']
            ])

    print(f"Analysis complete. Results written to {output_file}")
    print(f"Excluded proteins written to {excluded_file}")

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--run-analysis':
        main()
    else:
        python_path = create_virtual_environment()
        run_in_virtual_environment(python_path)