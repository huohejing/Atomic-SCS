#!/usr/bin/env python3
"""
Generate comprehensive validation sets:
- normal_100.smi: 100 molecules that are chemically valid (including small rings with strain)
- abnormal_100.smi: 100 molecules with various violations and edge cases
"""

import random

def generate_normal_molecules(n=100):
    """Normal molecules: chemically valid, may include small rings (cyclopropane, cyclobutane, etc.)"""
    base = [
        # Alkanes (linear & branched)
        "C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC",
        "CC(C)C", "CC(C)(C)C", "CC(C)CC", "CC(C)C(C)C", "CC(C)CC(C)C",
        # Cycloalkanes (including small rings)
        "C1CC1",        # cyclopropane (strain, but valid)
        "C1CCC1",       # cyclobutane (strain)
        "C1CCCC1",      # cyclopentane (minimal strain)
        "C1CCCCC1",     # cyclohexane (strain-free)
        "C1CCCCCC1",    # cycloheptane
        "C1CCCCCCC1",   # cyclooctane
        # Alkenes
        "C=C", "C=CC", "C=CCC", "C=CCCC", "CC=C(C)C", "CC=C", "C=C(C)C",
        # Alkynes
        "C#C", "C#CC", "C#CCC", "C#CCCC", "CC#C", "CC#CC",
        # Aromatics
        "c1ccccc1", "Cc1ccccc1", "c1ccc2ccccc2c1", "c1ccccc1C", "c1ccncc1",
        "c1cnccn1", "c1c[nH]cn1", "c1ccc2c(c1)cccc2", "c1ccc2ccccc2c1C",
        # Alcohols
        "CO", "CCO", "CC(C)O", "CCCCO", "C1CCCC1O", "c1ccccc1O", "CC(C)(C)O",
        # Ethers
        "COC", "CCOC", "CC(C)OC", "c1ccccc1OC", "CCOCC", "COCCOC",
        # Aldehydes
        "CC=O", "CCC=O", "CC(C)C=O", "c1ccccc1C=O",
        # Ketones
        "CC(=O)C", "CC(=O)CC", "CCC(=O)CC", "c1ccccc1C(=O)C",
        # Carboxylic acids
        "CC(=O)O", "CCC(=O)O", "CC(C)C(=O)O", "c1ccccc1C(=O)O",
        # Esters
        "CC(=O)OC", "CCOC(=O)C", "CC(=O)OCC", "c1ccccc1C(=O)OC",
        # Amines
        "CN", "CCN", "CC(C)N", "CNC", "CN(C)C", "c1ccccc1N", "CCNCC",
        # Amides
        "CC(=O)N", "CC(=O)NC", "CC(=O)N(C)C", "c1ccccc1C(=O)N",
        # Thiols & sulfides
        "CS", "CCS", "CSC", "CCSC", "CSCCS",
        # Heterocycles (including small rings)
        "C1OC1", "C1COC1", "C1CNC1", "C1CNCC1", "C1CCNCC1",
        "C1=CNC=C1", "C1=NC=CC=C1", "C1=CN=CC=C1", "C1=CN=CC=N1",
        # Additional drug-like fragments
        "CC(=O)Nc1ccccc1", "CC(=O)OCCN", "CN1CCCC1", "CC1CCCCC1",
        "C1CC2CCC1CC2", "C1CC2CC1CC2", "C1CC2CCCC12", "C1CC2CCC1C2",
    ]
    base = list(dict.fromkeys(base))
    random.shuffle(base)
    # Ensure exactly n
    while len(base) < n:
        base.append("C" * (len(base) % 10 + 1))
    return base[:n]

def generate_abnormal_molecules(n=100):
    """Problematic molecules: valence issues, extreme charges, invalid SMILES, etc."""
    problematic = []

    # 1. Valence violations
    valence = [
        "[C+](C)(C)(C)(C)", "[C+](C)(C)(C)(C)(C)", "[C+5]",
        "CC(C)(C)(C)(C)C", "C(C)(C)(C)(C)C",
        "[N+](C)(C)(C)(C)", "[N+](C)(C)(C)(C)(C)", "[N+5]",
        "[O+](C)(C)(C)", "[O+](C)(C)(C)(C)",
        "[P+](C)(C)(C)(C)", "[P+5]",
        "[S+](C)(C)(C)", "[S+6]",
    ]
    problematic.extend(valence[:15])

    # 2. Extreme charges
    charge = [
        "[C++]", "[C+++]", "[C--]", "[C---]",
        "[N++]", "[N--]",
        "[O++]", "[O--]",
        "[F++]", "[F--]",
        "[Cl++]", "[Cl--]",
        "[Li--]",
    ]
    problematic.extend(charge[:10])

    # 3. Small ring strain (already in normal? but we include extra examples)
    ring_strain = [
        "C1C1", "C1CC1", "C1CCC1", "C1OC1", "C1SC1", "C1NC1", "C1NCC1",
        "C1CC1C1", "C1CCCC1C1", "C1CCC1C1", "C1CC1CC1",
    ]
    problematic.extend(ring_strain[:10])

    # 4. Aromaticity issues
    aromatic = [
        "C1=CC=CC1",       # cyclopentadiene
        "C1=CC=C1",        # cyclobutadiene
        "C1=CCCC1",        # cyclopentene
        "C1=CCC=C1",       # 1,3-cyclohexadiene
        "c1c[nH+]cc1",     # pyridinium
    ]
    problematic.extend(aromatic[:5])

    # 5. Stereochemistry issues
    stereo = [
        "C[C@H](C)C", "C[C@](C)(C)C", "C[C@H](C)[C@H](C)C",
        "C[C@H](O)[C@H](O)C", "C[C@H](Cl)[C@@H](Cl)C",
    ]
    problematic.extend(stereo[:5])

    # 6. Edge cases (model limitations)
    edge = [
        "[BH4-]", "[BH3]C", "[B(OH)3]", "[SiH4]", "[SiH3]C",
        "[GeH4]", "[SeH2]", "[TeH2]", "[AsH3]", "[SbH3]",
        "CS(=O)C", "CS(=O)(=O)C", "CP(=O)(OC)OC", "COP(=O)(OC)OC",
        "[NH4+]", "[NH3+]", "[PH4+]", "[SH3+]", "[OH3+]",
        "[C+]", "[C+](C)(C)(C)C",
        "[NH3+][CH2]C(=O)[O-]",
        "C1CCCCCCCCCCCCCCCC1",
    ]
    problematic.extend(edge[:20])

    # 7. Invalid SMILES (parse errors)
    invalid = [
        "C(C", "C1CC1C1", "CC(C)(C)(C", "C1C", "[C+++]", "[C---]",
        "C=C=C=C", "C1CCC1C1", "C1C2C3C4", "C1C2C1C2",
        "invalid_smiles", "C1CC", "CC(C", "C1CC1C", "C1C1C1",
    ]
    problematic.extend(invalid[:15])

    # Remove duplicates and ensure exactly n
    problematic = list(dict.fromkeys(problematic))
    random.shuffle(problematic)
    if len(problematic) < n:
        filler = ["C1CC1", "[C+5]", "C(C)(C)(C)(C)C"]
        while len(problematic) < n:
            problematic.append(filler[len(problematic) % len(filler)])
    return problematic[:n]

def write_smiles(smiles_list, filename):
    with open(filename, 'w') as f:
        for s in smiles_list:
            f.write(s + "\n")

def main():
    normal = generate_normal_molecules(100)
    abnormal = generate_abnormal_molecules(100)
    write_smiles(normal, "normal_100.smi")
    write_smiles(abnormal, "abnormal_100.smi")
    print(f"Generated {len(normal)} normal and {len(abnormal)} abnormal SMILES.")
    print("Normal first 5:", normal[:5])
    print("Abnormal first 5:", abnormal[:5])

if __name__ == "__main__":
    main()