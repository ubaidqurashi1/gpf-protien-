from Bio.SeqUtils import GC
from Bio.Seq import Seq
import re

def dna_from_protein_seq(protein_seq, organism='ecoli'):
    """Convert protein to DNA using preferred codons (simplified)"""
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
