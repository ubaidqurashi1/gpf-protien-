#!/bin/bash

# Create GPF Protein Design Repository
# Author: You
# Date: 2025-04-05

set -e  # Exit on error

REPO_NAME="gpf-protein-design"
echo "🚀 Creating GPF repository: $REPO_NAME"

mkdir -p "$REPO_NAME"
cd "$REPO_NAME"

# --- LICENSE ---
cat > LICENSE << 'EOF'
MIT License

Copyright (c) 2025 Your Name

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
EOF

# --- requirements.txt ---
cat > requirements.txt << 'EOF'
biopython>=1.78
numpy>=1.21
scikit-learn>=1.0
matplotlib>=3.5
pandas>=1.3
pytest>=6.0
EOF

# --- gpf/ ---
mkdir -p gpf

cat > gpf/__init__.py << 'EOF'
from .core import gpf_transform_v4
from .predictors import predict_tm_thermonet_lite, predict_solubility, predict_expression
EOF

cat > gpf/utils.py << 'EOF'
from Bio.SeqUtils import GC
from Bio.Seq import Seq
import re

def dna_from_protein_seq(protein_seq, organism='ecoli'):
    ecoli_codons = {
        'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
        'Q': 'CAG', 'E': 'GAA', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
        'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTC', 'P': 'CCG',
        'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTG'
    }
    dna = ''.join(ecoli_codons.get(aa, 'NNN') for aa in protein_seq)
    return dna

def extract_expression_features(promoter_seq=None, utr5_seq=None, utr3_seq=None):
    features = {}
    features['promoter_strength'] = promoter_seq.count('TATAAT') * 2 + promoter_seq.count('TTGACA') * 2 if promoter_seq else 0
    features['rbs_score'] = len(re.findall(r'[AG]G[AG]G', utr5_seq[-20:])) if utr5_seq else 0
    features['utr5_gc'] = GC(Seq(utr5_seq)) if utr5_seq else 50
    features['utr5_length'] = len(utr5_seq) if utr5_seq else 0
    features['polyA_signal'] = utr3_seq.count('AATAAA') if utr3_seq else 0
    features['utr3_length'] = len(utr3_seq) if utr3_seq else 0
    return [features['promoter_strength'], features['rbs_score'], features['utr5_gc'],
            features['utr5_length'], features['polyA_signal'], features['utr3_length']]
EOF

cat > gpf/predictors.py << 'EOF'
import numpy as np

def predict_tm_thermonet_lite(protein_seq):
    if not protein_seq:
        return 50.0
    L = len(protein_seq)
    kd = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
          'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
    helix_prop = {'A':1.45,'R':1.00,'N':0.76,'D':0.92,'C':0.77,'Q':1.17,'E':1.53,'G':0.53,'H':1.24,
                  'I':1.00,'L':1.34,'K':1.07,'M':1.20,'F':1.12,'P':0.59,'S':0.79,'T':0.82,'W':1.14,'Y':0.61,'V':1.00}
    hydro = np.mean([kd.get(aa, 0) for aa in protein_seq])
    helix = np.mean([helix_prop.get(aa, 1.0) for aa in protein_seq])
    arom = (protein_seq.count('F') + protein_seq.count('W') + protein_seq.count('Y')) / L
    charged = sum(protein_seq.count(aa) for aa in 'RKDE') / L
    pro_loop = protein_seq.count('P') * 0.5 / L
    cys = protein_seq.count('C') / L
    agg = sum(protein_seq.count(aa) for aa in 'VILFWY') / L
    tm_pred = (45.0 + 3.2 * hydro + 12.0 * (helix - 1.0) + 8.0 * arom +
               5.0 * charged + 4.0 * pro_loop + 6.0 * cys - 15.0 * agg)
    return np.clip(tm_pred, 30, 95)

def predict_solubility(protein_seq):
    kd = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
          'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
    hydro = np.mean([kd.get(aa, 0) for aa in protein_seq])
    return 1.0 / (1.0 + np.exp(2 * (hydro - 0)))

def predict_expression(expr_features):
    return np.clip(0.2 * expr_features[0] + 0.5 * expr_features[1] + 0.1 * (100 - expr_features[2]), 0, 10)
EOF

cat > gpf/core.py << 'EOF'
import numpy as np
from .utils import extract_expression_features
from .predictors import predict_tm_thermonet_lite, predict_solubility, predict_expression

def gpf_transform_v4(protein_seq, promoter_seq=None, utr5_seq=None, utr3_seq=None, experimental_data=None):
    if not protein_seq:
        raise ValueError("Protein sequence cannot be empty")
    
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    aa_count = np.array([protein_seq.count(aa) for aa in aa_list])
    aa_freq = aa_count / len(protein_seq)
    
    kd = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
          'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
    hydro = np.mean([kd.get(aa, 0) for aa in protein_seq])
    basic = sum(protein_seq.count(aa) for aa in 'RKH')
    acidic = sum(protein_seq.count(aa) for aa in 'DE')
    net_charge = basic - acidic
    aromatic_frac = sum(protein_seq.count(aa) for aa in 'FWY') / len(protein_seq)
    disorder_score = sum(protein_seq.count(aa) for aa in 'PESQ') / len(protein_seq)
    length = len(protein_seq)
    helix_prop = {'A':1.45,'R':1.00,'N':0.76,'D':0.92,'C':0.77,'Q':1.17,'E':1.53,'G':0.53,'H':1.24,
                  'I':1.00,'L':1.34,'K':1.07,'M':1.20,'F':1.12,'P':0.59,'S':0.79,'T':0.82,'W':1.14,'Y':0.61,'V':1.00}
    helix = np.mean([helix_prop.get(aa, 1.0) for aa in protein_seq])
    strand = 0.3
    coil = 1 - helix - strand
    rsa = sum(protein_seq.count(aa) for aa in 'RNDQEHKSTY') / len(protein_seq)
    agg = sum(protein_seq.count(aa) for aa in 'VILFWY') / len(protein_seq)
    protein_feat = np.array([hydro, net_charge, aromatic_frac, disorder_score, length, helix, strand, coil, rsa, agg])
    
    expr_feat = extract_expression_features(promoter_seq, utr5_seq, utr3_seq)
    
    prior_tm = predict_tm_thermonet_lite(protein_seq)
    prior_sol = predict_solubility(protein_seq)
    prior_expr = predict_expression(expr_feat)
    
    if experimental_
        obs_tm = experimental_data.get('Tm', prior_tm)
        obs_sol = experimental_data.get('solubility', prior_sol)
        obs_expr = experimental_data.get('expression', prior_expr)
        posterior_tm = (prior_tm + obs_tm) / 2
        posterior_sol = 0.6 * prior_sol + 0.4 * obs_sol
        posterior_expr = 0.7 * prior_expr + 0.3 * obs_expr
    else:
        posterior_tm, posterior_sol, posterior_expr = prior_tm, prior_sol, prior_expr
    
    Z = np.concatenate([aa_freq, protein_feat, expr_feat, [posterior_tm, posterior_sol, posterior_expr]])
    return Z
EOF

# --- data/ ---
mkdir -p data

cat > data/gfp_mutants.csv << 'EOF'
name,position,wt_aa,mut_aa,sequence
WT,,,,"MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
V163D,163,V,D,"MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLDNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
L201K,201,L,K,"MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIEKKIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
F165E,165,F,E,"MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEEVTAAGITLGMDELYK"
I180R,180,I,R,"MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGRTLGMDELYK"
V217D,217,V,D,"MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
EOF

# --- notebooks/ ---
mkdir -p notebooks
echo "# GPF Demo Notebook" > notebooks/01_gpf_demo.ipynb
echo "# Analysis Notebook" > notebooks/02_analysis.ipynb

# --- results/ ---
mkdir -p results
touch results/.gitkeep

# --- experiments/ ---
mkdir -p experiments

cat > experiments/protocol.md << 'EOF'
# Experimental Validation Protocol: GPF-Designed GFP Mutants

## Cloning
- **Vector**: pET28a(+)
- **Insert**: GFP variants (WT + 5 mutants)
- **Method**: Gibson Assembly or NdeI/XhoI digest
- **Verification**: Sanger sequencing

## Expression
- **Strain**: E. coli BL21(DE3)
- **Culture**: LB + 50 µg/mL kanamycin
- **Induction**: OD600 = 0.6 → 0.5 mM IPTG, 18°C, 16h

## Assays
1. **Fluorescence**: 
   - Plate reader: Ex 395 nm / Em 509 nm
   - Normalize to OD600

2. **Thermal Stability (DSF)**:
   - Lysate + 5X SYPRO Orange
   - Ramp: 25–95°C at 1°C/min
   - Tm = inflection point

3. **Solubility**:
   - Lyse 1 mL culture
   - Centrifuge 15,000g, 10 min
   - Run soluble (supernatant) and insoluble (pellet) on 12% SDS-PAGE
   - Quantify band intensity (ImageJ)

## Success Criteria
- ΔTm ≥ –3.0°C vs WT
- Solubility ≥ 85%
- Fluorescence ≥ 90% of WT
EOF

# --- tests/ (UNIT TESTS) ---
mkdir -p tests

cat > tests/test_gpf.py << 'EOF'
import pytest
import numpy as np
from gpf import gpf_transform_v4

def test_gpf_wt_gfp():
    wt_seq = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    Z = gpf_transform_v4(wt_seq)
    assert len(Z) == 39
    assert 60 < Z[-3] < 70  # Tm between 60-70°C
    assert 0.6 < Z[-2] < 0.9  # Solubility reasonable

def test_gpf_v163d():
    wt_seq = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    mut_seq = wt_seq[:162] + "D" + wt_seq[163:]
    Z_wt = gpf_transform_v4(wt_seq)
    Z_mut = gpf_transform_v4(mut_seq)
    assert Z_mut[-2] > Z_wt[-2]  # Solubility should increase
    assert abs(Z_mut[-3] - Z_wt[-3]) < 5  # Tm change mild

def test_empty_sequence():
    with pytest.raises(ValueError):
        gpf_transform_v4("")
EOF

# --- tools/foldx/ (FOLDX INTEGRATION) ---
mkdir -p tools/foldx

cat > tools/foldx/README.md << 'EOF'
# FoldX Integration for GPF

## Requirements
- FoldX v5+ installed and in PATH
- PDB file (e.g., 1GFL.pdb)

## Usage
```bash
./run_foldx.sh 1GFL.pdb V163D
