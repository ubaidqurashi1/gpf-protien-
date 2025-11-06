# gpf/core.py
import numpy as np
from .utils import extract_expression_features

def gpf_transform_v5(protein_seq, promoter_seq=None, utr5_seq=None, utr3_seq=None, experimental_data=None):
    """
    GPF v5: Fully aligned with final manuscript.
    """
    if not protein_seq:
        raise ValueError("Protein sequence cannot be empty")
    
    L = len(protein_seq)
    
    # === 1. Amino Acid Frequency (20D) ===
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    aa_count = np.array([protein_seq.count(aa) for aa in aa_list])
    f_aa = aa_count / L
    
    # === 2. Protein Features (10D) ===
    # Hydrophobicity (Kyte-Doolittle)
    kd_scale = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
                'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
    H = np.mean([kd_scale.get(aa, 0) for aa in protein_seq])
    
    # Net charge
    Q = (protein_seq.count('K') + protein_seq.count('R') + protein_seq.count('H') - 
         protein_seq.count('D') - protein_seq.count('E'))
    
    # Aromatic fraction
    A = (protein_seq.count('F') + protein_seq.count('W') + protein_seq.count('Y')) / L
    
    # Disorder propensity
    D = (protein_seq.count('P') + protein_seq.count('E') + protein_seq.count('S') + protein_seq.count('Q')) / L
    
    # Length
    length = L
    
    # Secondary structure (Chou-Fasman)
    helix_prop = {'A':1.45,'R':1.00,'N':0.76,'D':0.92,'C':0.77,'Q':1.17,'E':1.53,'G':0.53,'H':1.24,
                  'I':1.00,'L':1.34,'K':1.07,'M':1.20,'F':1.12,'P':0.59,'S':0.79,'T':0.82,'W':1.14,'Y':0.61,'V':1.00}
    strand_prop = {'A':0.99,'R':0.93,'N':0.76,'D':0.61,'C':1.30,'Q':1.43,'E':1.67,'G':0.83,'H':1.67,
                   'I':1.60,'L':1.22,'K':1.15,'M':1.45,'F':1.90,'P':0.55,'S':0.83,'T':1.20,'W':1.90,'Y':1.70,'V':1.52}
    alpha = np.mean([helix_prop.get(aa, 1.0) for aa in protein_seq])
    beta = np.mean([strand_prop.get(aa, 1.0) for aa in protein_seq])
    gamma = 1.0 - alpha - beta
    
    # Solvent accessibility proxy
    polar_residues = "RNDQEHKSTY"
    R = sum(protein_seq.count(aa) for aa in polar_residues) / L
    
    # Aggregation propensity
    agg_residues = "VILFWY"
    G = sum(protein_seq.count(aa) for aa in agg_residues) / L
    
    f_prot = np.array([H, Q, A, D, length, alpha, beta, gamma, R, G])
    
    # === 3. Regulatory Features (6D) ===
    f_expr = np.array(extract_expression_features(promoter_seq, utr5_seq, utr3_seq))
    
    # === 4. Priors ===
    # Stability prior (Equation 4)
    charged_count = sum(protein_seq.count(aa) for aa in 'RKDEH')
    pro_count = protein_seq.count('P')
    cys_count = protein_seq.count('C')
    Tm_prior = (45.0 + 
                3.2 * H + 
                12.0 * (alpha - 1.0) + 
                8.0 * A + 
                5.0 * (charged_count / L) + 
                4.0 * (pro_count / (2 * L)) + 
                6.0 * (cys_count / L) - 
                15.0 * G)
    Tm_prior = np.clip(Tm_prior, 30, 95)
    
    # Solubility prior (Equation 5)
    S_prior = 1.0 / (1.0 + np.exp(2.0 * (H - 0.0)))
    
    # Expression prior (Equation 6)
    P, RBS, utr_gc, utr5_len, polyA, utr3_len = f_expr
    E_prior = 0.2 * P + 0.5 * RBS + 0.1 * (100 - utr_gc)
    E_prior = min(10.0, E_prior)
    
    # === 5. Bayesian Updating ===
    if experimental_data:
        Tm_obs = experimental_data.get('Tm', Tm_prior)
        S_obs = experimental_data.get('solubility', S_prior)
        E_obs = experimental_data.get('expression', E_prior)
        
        # Uncertainties
        sigma_T_prior, sigma_T_obs = 5.0, 3.0
        sigma_S_prior, sigma_S_obs = 0.15, 0.10
        sigma_E_prior, sigma_E_obs = 2.0, 1.5
        
        Tm_post = (Tm_prior / sigma_T_prior**2 + Tm_obs / sigma_T_obs**2) / (1/sigma_T_prior**2 + 1/sigma_T_obs**2)
        S_post = (S_prior / sigma_S_prior**2 + S_obs / sigma_S_obs**2) / (1/sigma_S_prior**2 + 1/sigma_S_obs**2)
        E_post = (E_prior / sigma_E_prior**2 + E_obs / sigma_E_obs**2) / (1/sigma_E_prior**2 + 1/sigma_E_obs**2)
    else:
        Tm_post, S_post, E_post = Tm_prior, S_prior, E_prior
    
    # === 6. Full 39D Percept ===
    Z = np.concatenate([f_aa, f_prot, f_expr, [Tm_post, S_post, E_post]])
    assert len(Z) == 39, f"Z has {len(Z)} dimensions, expected 39"
    
    return Z
