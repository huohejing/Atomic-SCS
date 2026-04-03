"""
Molecule Analyzer Module
Analyzes molecular structure and computes compliance scores
"""

from rdkit import Chem
from atomic_rules import get_all_rules


class MoleculeAnalyzer:
    def __init__(self):
        self.rules = get_all_rules()
    
    def analyze_molecule(self, mol, strictness='balanced'):
        """
        Analyze a molecule and return detailed compliance analysis.
        
        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit molecule object.
        strictness : str, optional
            'conservative', 'balanced', or 'liberal'. Default 'balanced'.
        
        Returns
        -------
        dict
            Analysis results including scores for each atom.
        """
        if mol is None:
            return {
                'valid': False,
                'error': 'Invalid molecule',
                'atom_details': [],
                'num_atoms': 0
            }
        
        # Compute ring strain scores at molecular level based on strictness
        ring_strain_scores = self._compute_ring_strain_scores(mol, strictness)
        
        atom_analysis = []
        for i, atom in enumerate(mol.GetAtoms()):
            atom_data = {
                'index': i,
                'element': atom.GetSymbol(),
                'formal_charge': atom.GetFormalCharge(),
                'degree': atom.GetDegree(),
                'hybridization': str(atom.GetHybridization()),
                'is_aromatic': atom.GetIsAromatic(),
                'scores': {},
                'violations': []
            }
            
            # Apply all rules to this atom
            for rule_name, rule_func in self.rules.items():
                if rule_name == 'ring_strain_check':
                    # Use precomputed ring strain score instead of calling function
                    score = ring_strain_scores.get(i, 0.0)
                    is_valid = (score == 0.0)
                    explanation = f"Ring strain score: {score:.3f}"
                else:
                    # Call rule function with strictness parameter
                    is_valid, score, explanation = rule_func(atom, strictness)
                
                atom_data['scores'][rule_name] = score
                if not is_valid:
                    atom_data['violations'].append({
                        'rule': rule_name,
                        'score': score,
                        'explanation': explanation
                    })
            
            atom_analysis.append(atom_data)
        
        return {
            'valid': True,
            'atom_details': atom_analysis,
            'num_atoms': mol.GetNumAtoms()
        }
    
    def _compute_ring_strain_scores(self, mol, strictness):
        """
        Compute ring strain scores at molecular level to avoid double counting.
        Each ring is scored once and the score is distributed among ring atoms.
        Strictness affects the penalty strength.
        
        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit molecule object.
        strictness : str
            'conservative', 'balanced', or 'liberal'.
        
        Returns
        -------
        dict
            {atom_index: ring_strain_score}
        """
        ring_strain_scores = {}
        
        # Define penalty strengths based on strictness
        ring_penalties = {
            'conservative': {
                3: 0.3,  # Cyclopropane - higher penalty in conservative mode
                4: 0.15  # Cyclobutane - higher penalty in conservative mode
            },
            'balanced': {
                3: 0.2,  # Cyclopropane - standard penalty
                4: 0.1   # Cyclobutane - standard penalty
            },
            'liberal': {
                3: 0.1,  # Cyclopropane - lower penalty in liberal mode
                4: 0.05  # Cyclobutane - lower penalty in liberal mode
            }
        }
        
        # Get penalties for current strictness
        penalties = ring_penalties.get(strictness, ring_penalties['balanced'])
        
        # Get all simple rings
        sssr = mol.GetRingInfo().AtomRings()
        
        # Track which rings each atom belongs to
        atom_to_rings = {}
        for ring_idx, ring_atoms in enumerate(sssr):
            ring_size = len(ring_atoms)
            
            if ring_size in penalties:
                ring_penalty = penalties[ring_size]
                penalty_per_atom = ring_penalty / ring_size
                
                for atom_idx in ring_atoms:
                    if atom_idx not in atom_to_rings:
                        atom_to_rings[atom_idx] = []
                    atom_to_rings[atom_idx].append((ring_idx, penalty_per_atom))
        
        # Assign scores ensuring no double counting
        for atom_idx, rings_info in atom_to_rings.items():
            # For atoms in multiple rings, take the maximum penalty rather than summing
            # This prevents over-penalizing bridged ring systems
            max_penalty = max(penalty for _, penalty in rings_info)
            ring_strain_scores[atom_idx] = max_penalty
        
        return ring_strain_scores


def analyze_single_molecule(smiles, strictness='balanced'):
    """
    Convenience function to analyze a single molecule from SMILES string.
    
    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    strictness : str, optional
        'conservative', 'balanced', or 'liberal'. Default 'balanced'.
    
    Returns
    -------
    dict
        Analysis results.
    """
    mol = Chem.MolFromSmiles(smiles)
    analyzer = MoleculeAnalyzer()
    return analyzer.analyze_molecule(mol, strictness)


if __name__ == "__main__":
    # Test the analyzer
    test_smiles = [
        "C",           # Methane
        "C1CC1",       # Cyclopropane
        "C1CCC1",      # Cyclobutane
        "C1CCCC1",     # Cyclopentane
        "C1CCCCCC1",   # Cycloheptane
        "[NH4+]",      # Ammonium
        "[C+]",        # Carbocation
        "CC(=O)[O-]",  # Acetate
    ]
    
    analyzer = MoleculeAnalyzer()
    
    for smiles in test_smiles:
        print(f"\nAnalyzing: {smiles}")
        result = analyzer.analyze_molecule(Chem.MolFromSmiles(smiles))
        print(f"Number of atoms: {result['num_atoms']}")
        if result['valid']:
            for atom_data in result['atom_details']:
                if atom_data['violations']:
                    print(f"  Atom {atom_data['index']} ({atom_data['element']}): {atom_data['violations']}")