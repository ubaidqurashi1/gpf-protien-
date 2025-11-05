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
