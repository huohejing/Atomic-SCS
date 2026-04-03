# Atomic-SCS: Atomic‑level Chemical Rule Scoring

**Atomic-SCS** is an open‑source Python tool for atomic‑level chemical rule evaluation. It assesses each atom in a molecule across multiple dimensions (valence, charge, aromaticity, ring strain, ring count) and returns a continuous compliance score (0 = fully compliant, 1 = severe violation). The tool provides three strictness levels, three output modes (assess, diagnose, repair), and supports multiple output formats (text/TSV, CSV, JSON). Its continuous scores are differentiable, making it suitable as a constraint or reward signal in generative molecular design pipelines.

---

## ✨ Features

- **Atom‑level multi‑dimensional scoring**  
  Valence (quadratic penalty), charge (distance penalty), aromaticity, ring strain (3‑/4‑membered rings), ring count.
- **Three strictness modes**  
  `conservative`, `balanced`, `liberal` – adjust compliance thresholds and hypervalence tolerance.
- **Three operation modes**  
  - `assess` – summary of molecular score, compliance, confidence level.  
  - `diagnose` – detailed per‑atom scores and issues.  
  - `repair` – prioritized fix suggestions.
- **Multiple output formats**  
  `text` (default TSV), `csv`, `json` (including detailed per‑atom verbose output).
- **Command‑line interface** and **Python API**.
- **Comprehensive validation** – 100 normal vs. 100 problematic molecules, Mann‑Whitney U test p < 1e-18.

---

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
Argument	Choices	Default	Description
input_file	–	–	SMILES file (one per line)
--mode, -m	assess, diagnose, repair	assess	Operation mode
--strictness, -s	conservative, balanced, liberal	balanced	Strictness level
--verbose, -v	–	False	Detailed atom‑level output
--format, -f	text, csv, json	text	Output format
--output, -o	–	stdout	Output file
Output formats
text (TSV) – Tab‑separated values, easy to view or import.

csv – Comma‑separated, suitable for spreadsheets.

json – Structured data, ideal for programmatic processing.

🧪 Validation
We validated Atomic-SCS on 100 normal molecules (covering diverse chemical classes) and 100 problematic molecules (containing valence errors, extreme charges, ring strain, invalid SMILES, etc.). The normal group scores are tightly centered at 0 (with a few small‑ring outliers receiving minor ring‑strain penalties), while the problematic group shows a wide distribution. A Mann‑Whitney U test yields p < 1e-18, demonstrating excellent discriminative power.

Note: You can generate the boxplot by running python plot_from_tsv.py after scoring your validation sets.

📖 Citation
If you use Atomic-SCS in your research, please cite:

text
Author(s). (2026). Atomic-SCS: Atomic‑level Differentiable Chemical Rule Scoring for Molecular Design.
Journal / preprint, DOI: ... (to be added).
🤝 Contributing
Issues and pull requests are welcome. Please open an issue for discussion before major changes.

📄 License
Atomic-SCS is released under the MIT License. See LICENSE for details.

🙏 Acknowledgements
RDKit – cheminformatics toolkit.

Chemical rules based on standard textbooks: Clayden Organic Chemistry, IUPAC Blue Book, 邢其毅《基础有机化学》。
