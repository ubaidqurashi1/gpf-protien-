import numpy as np
from .utils import dna_from_protein_seq, extract_expression_features
from .predictors import predict_tm_thermonet_lite, predict_solubility, predict_expression

def gpf_transform_v4(protein_seq, promoter_seq=None, utr5_seq=None, utr3_seq=None, experimental_data=None):
    # Amino acid freq (20D)
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    aa_count = np.array([protein_seq.count(aa) for aa in aa_list])
    aa_freq = aa_count / (len(protein_seq) + 1e-8)
    
    # Protein features (10D)
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
    strand = 0.3  # placeholder
    coil = 1 - helix - strand
    rsa = sum(protein_seq.count(aa) for aa in 'RNDQEHKSTY') / len(protein_seq)
    agg = sum(protein_seq.count(aa) for aa in 'VILFWY') / len(protein_seq)
    protein_feat = np.array([hydro, net_charge, aromatic_frac, disorder_score, length, helix, strand, coil, rsa, agg])
    
    # Expression features (6D)
    expr_feat = extract_expression_features(promoter_seq, utr5_seq, utr3_seq)
    
    # Priors
    prior_tm = predict_tm_thermonet_lite(protein_seq)
    prior_sol = predict_solubility(protein_seq)
    prior_expr = predict_expression(expr_feat)
    
    # Bayesian update
    if experimental_data:
        obs_tm = experimental_data.get('Tm', prior_tm)
        obs_sol = experimental_data.get('solubility', prior_sol)
        obs_expr = experimental_data.get('expression', prior_expr)
        posterior_tm = (prior_tm + obs_tm) / 2
        posterior_sol = 0.6 * prior_sol + 0.4 * obs_sol
        posterior_expr = 0.7 * prior_expr + 0.3 * obs_expr
    else:
        posterior_tm, posterior_sol, posterior_expr = prior_tm, prior_sol, prior_expr
    
    # Full Z (20 + 10 + 6 + 3 = 39D)
    Z = np.concatenate([aa_freq, protein_feat, expr_feat, [posterior_tm, posterior_sol, posterior_expr]])
    return Z
