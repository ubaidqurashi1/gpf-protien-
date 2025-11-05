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
