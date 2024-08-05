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

def check_n_glycosylation_sites(sequence):
    pattern = re.compile(r'N[^P][ST]')
    return len(pattern.findall(sequence))

def analyze_protein(sequence, features, max_size=1500):
    if len(sequence) > max_size:
        return None, f"Protein too large ({len(sequence)} aa)", "Low", None

    compatible_ptms = {
        "phospho": "Phosphorylation",
        "n-linked glycosyl": "N-linked glycosylation",
        "glcnac-asparagine": "N-linked glycosylation",
        "oligosaccharide-asparagine": "N-linked glycosylation",
        "simple mannose-type": "Simple mannose-type glycan",
        "s-layer glyco": "S-layer glycoprotein",
        "o-linked glycosyl": "Simple O-linked glycosylation",
        "methyl": "Methylation",
        "acetyl": "Acetylation",
        "disulfide": "Disulfide bond"
    }

    excluded_ptms = [
        "complex glycosyl", "complex n-linked", "gal-galnac",
        "sialic acid", "fucosyl", "ubiquitin", "sumo", "farnesyl", 
        "geranylgeranyl", "palmitoyl", "myristoyl"
    ]

    found_compatible_ptms = set()
    n_glycosylation_sites = check_n_glycosylation_sites(sequence)
    
    ptms = [feature for feature in features if feature.type == "Site"]
    has_ptm_annotations = bool(ptms)

    for ptm in ptms:
        ptm_description = ptm.qualifiers.get('note', [''])[0].lower()
        if any(re.search(rf'\b{excluded_ptm}\w*', ptm_description) for excluded_ptm in excluded_ptms):
            if "not " in ptm_description or "absence of " in ptm_description:
                continue
            return None, f"Excluded PTM found: {ptm_description}", "Low", None
        
        for ptm_key, ptm_name in compatible_ptms.items():
            if ptm_key in ptm_description:
                found_compatible_ptms.add(ptm_name)
    
    cys_count = sequence.count('C')
    cys_density = (cys_count / len(sequence)) * 100
    
    pass_reasons = []
    if found_compatible_ptms:
        pass_reasons.append(f"Compatible PTMs found: {', '.join(found_compatible_ptms)}")
    elif has_ptm_annotations:
        pass_reasons.append("No incompatible PTMs found")
    else:
        pass_reasons.append("No PTM annotations available")
    
    if n_glycosylation_sites > 0:
        pass_reasons.append(f"Potential N-glycosylation sites: {n_glycosylation_sites}")
    
    if cys_density > 5:
        pass_reasons.append("High cysteine density")
    
    pass_reason = "; ".join(pass_reasons)
    
    confidence = "High" if has_ptm_annotations else "Medium"
    
    return {
        'length': len(sequence),
        'cys_density': cys_density,
        'n_glycosylation_sites': n_glycosylation_sites
    }, "Passed", confidence, pass_reason

def rank_proteins(accessions):
    results = []
    excluded = []
    for accession in accessions:
        try:
            name, sequence, features = fetch_protein_info(accession)
            analysis, message, confidence, pass_reason = analyze_protein(sequence, features)
            if analysis:
                analysis['accession'] = accession
                analysis['name'] = name
                analysis['screening_result'] = message
                analysis['confidence'] = confidence
                analysis['pass_reason'] = pass_reason
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

    return sorted(results, key=lambda x: x['cys_density'], reverse=True), excluded

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
        writer.writerow(['Accession', 'Name', 'Length', 'Cys Density', 'N-Glycosylation Sites', 'Screening Result', 'Confidence', 'Pass Reason'])
        for protein in ranked_proteins:
            writer.writerow([
                protein['accession'],
                protein['name'],
                protein['length'],
                f"{protein['cys_density']:.2f}",
                protein['n_glycosylation_sites'],
                protein['screening_result'],
                protein['confidence'],
                protein['pass_reason']
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
