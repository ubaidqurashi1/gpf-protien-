# recompute_predictions.py
import numpy as np
from gpf.core import gpf_transform_v5

# Wild-type GFP sequence (UniProt P42212)
WT_SEQ = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"

# Define mutants: (name, position_1_indexed, mutant_aa)
mutants = [
    ("V163D", 163, "D"),
    ("V217D", 217, "D"),
    ("L201K", 201, "K"),
    ("F165E", 165, "E"),
    ("I180R", 180, "R")
]

# Compute WT baseline
Z_wt = gpf_transform_v5(WT_SEQ)
Tm_wt = Z_wt[-3]
S_wt = Z_wt[-2]

print("GPF v5 Predictions (Updated Stability Model)")
print("=" * 50)
print(f"{'Mutant':<8} {'ΔTm (°C)':<10} {'ΔSolubility':<12} {'Tm (°C)':<8} {'Solubility':<10}")
print("-" * 50)

results = []
for name, pos, aa_mut in mutants:
    # Create mutant sequence (0-indexed)
    seq_list = list(WT_SEQ)
    seq_list[pos - 1] = aa_mut  # convert to 0-index
    mut_seq = "".join(seq_list)
    
    # Predict
    Z_mut = gpf_transform_v5(mut_seq)
    Tm_mut = Z_mut[-3]
    S_mut = Z_mut[-2]
    
    dTm = Tm_mut - Tm_wt
    dS = S_mut - S_wt
    
    print(f"{name:<8} {dTm:<10.1f} {dS:<12.3f} {Tm_mut:<8.1f} {S_mut:<10.3f}")
    results.append((name, dTm, dS, Tm_mut, S_mut))

# Save to CSV
import csv
with open('results/gpf_v5_predictions.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['Mutant', 'Delta_Tm_C', 'Delta_Solubility', 'Tm_C', 'Solubility'])
    for row in results:
        writer.writerow(row)

print("\nResults saved to results/gpf_v5_predictions.csv")
