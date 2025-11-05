# Genomic Perception Fusion (GPF)

A lightweight, interpretable framework for protein functional tuning—from DNA sequence to high-order functional representation.

## Features
- Predicts **stability (Tm)**, **solubility**, and **expression** from sequence
- Models **non-coding regions** (promoters, UTRs)
- **No GPU required** — runs on CPU in <1 sec
- Validated on **5 GFP surface mutants**

## Quick Start
```bash
pip install -r requirements.txt
python -c "from gpf import gpf_transform_v4; print(gpf_transform_v4('MASKGE...'))"
