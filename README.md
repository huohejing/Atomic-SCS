# Atomic-SCS: Atomic‑level Chemical Rule Scoring for Generative Molecular Design

Atomic-SCS is an open‑source Python tool that evaluates the chemical validity of molecules at the atom level. It scores each atom across four dimensions (valence, charge, aromaticity, ring strain) and returns a continuous compliance score (0 = fully compliant, 1 = severe violation). The tool supports three strictness levels (conservative, balanced, liberal) and three operation modes (assess, diagnose, repair). Output can be generated as text (TSV), CSV, or JSON. The continuous scores can serve as a reward signal in generative model training.

## ✨ Features

**Four scoring dimensions**

- Valence (quadratic penalty)
- Charge (distance penalty)
- Aromaticity (RDKit flags)
- Ring strain (only 3‑ and 4‑membered rings)

**Three strictness modes**  
`conservative`, `balanced`, `liberal` – adjust compliance thresholds and hypervalence tolerance for S and P.

**Three operation modes**

- `assess` – molecular score, compliance status, confidence level.
- `diagnose` – per‑atom scores and issue descriptions.
- `repair` – prioritized text suggestions (no automatic structure modification).

**Output formats**

- text (default TSV)
- csv
- json (including verbose per‑atom output)

**Command‑line interface and Python API**

## 🔧 Installation

### Using conda (recommended)

```bash
# Clone the repository
git clone https://github.com/huohejing/Atomic-SCS.git
cd Atomic-SCS

# Create and activate the conda environment
conda env create -f env.yml
conda activate atomic-scs
Manual installation
bash
conda create -n atomic-scs python=3.10 -y
conda activate atomic-scs
conda install -c conda-forge rdkit=2023.09.4 numpy pillow pandas matplotlib seaborn scipy -y
🚀 Quick Start
Command line
bash
# Assess a file of SMILES (TSV output)
python atomic_scs.py input.smi --mode assess --strictness balanced

# Detailed diagnosis (verbose)
python atomic_scs.py input.smi --mode diagnose --verbose

# Repair guidance (summary)
python atomic_scs.py input.smi --mode repair

# CSV output
python atomic_scs.py input.smi --mode assess --format csv

# JSON output
python atomic_scs.py input.smi --mode assess --format json
Python API
python
from atomic_scs import AtomicSCS

scs = AtomicSCS(strictness='balanced')
result = scs.assess_molecule("C1CC1")   # cyclopropane
print(f"Score: {result['score']:.3f}, compliant: {result['compliant']}")

# Detailed diagnosis
diagnosis = scs.diagnose_molecule("C1CC1")
print(diagnosis)

# Repair guidance
repairs = scs.repair_molecule("C1CC1")
print(repairs)
📚 Usage
Command‑line arguments

input_file : SMILES file (one per line)

--mode, -m : assess, diagnose, repair (default: assess) – Operation mode

--strictness, -s : conservative, balanced, liberal (default: balanced) – Strictness level

--verbose, -v : Detailed atom‑level output (default: False)

--format, -f : text, csv, json (default: text) – Output format

--output, -o : Output file (default: stdout)

Output formats

text (TSV) – Tab‑separated values, easy to view or import.

csv – Comma‑separated, suitable for spreadsheets.

json – Structured data, ideal for programmatic processing.

🧪 Validation and Performance
We validated Atomic-SCS on 100 normal molecules (covering diverse chemical classes) and 100 problematic molecules (containing valence errors, extreme charges, ring strain, invalid SMILES, etc.). The normal group scores are tightly centered at 0 (with small‑ring outliers receiving minor ring‑strain penalties), while the problematic group shows a wide distribution. A Mann‑Whitney U test yields p < 1e‑18, demonstrating excellent discriminative power.

Performance: On a single core of an Intel Core i7‑14650HX CPU (Ubuntu 22.04 via WSL), scoring 100 molecules (average ~30 heavy atoms) took 0.0071 seconds total (0.000071 seconds per molecule). The tool scales easily via multiprocessing.

⚠️ Limitations
Only common elements (H, C, N, O, F, P, S, Cl, Br, I) are fully supported; other elements default to a maximum valence of 4.

Aromaticity relies on RDKit's built‑in algorithm, which may be inaccurate for some fused ring systems.

Ring strain is only penalized for 3‑ and 4‑membered rings.

Valence check uses atom.GetDegree() (number of bonded heavy atoms) rather than total valence including hydrogens. This simplification is suitable for typical organic molecules; for molecules with unusual hydrogen counts, scoring may be less accurate.

Repair mode outputs only text suggestions and does not automatically modify molecular structures.

Installation must be done from source via a conda environment; pip installation is not yet supported.

Stereochemistry (R/S, E/Z) is ignored; a molecule with wrong chirality will not be penalized.

Invalid SMILES syntax (e.g., unmatched parentheses) causes immediate failure, returning a score of 1.0.

Resonance effects are not modeled; carboxylate anions and similar species may receive minor penalties due to formal charge deviations.

📖 Citation
If you use Atomic-SCS in your research, please cite the preprint:

Hejing Huo, Miaomiao Niu. Atomic-SCS: An Atom‑Level Chemical Rule Scoring Tool for Generative Molecular Design. Preprints.org 2026, 202604.0809.v1. DOI: 10.20944/preprints202604.0809.v1

🤝 Contributing
Issues and pull requests are welcome. Please open an issue for discussion before major changes.

📄 License
Atomic-SCS is released under the MIT License. See LICENSE for details.

🙏 Acknowledgements
RDKit – cheminformatics toolkit.

Chemical rules based on standard textbooks: Clayden Organic Chemistry, IUPAC Blue Book, 邢其毅《基础有机化学》。
