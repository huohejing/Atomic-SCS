"""
Atomic Rules Module
Defines detailed checks for individual atomic properties
"""

def check_valency_detailed(atom, strictness='balanced'):
    """
    Check if atom's degree matches its expected valency.
    
    Parameters
    ----------
    atom : rdkit.Chem.rdchem.Atom
        RDKit atom object.
    strictness : str, optional
        'conservative', 'balanced', or 'liberal'. Default 'balanced'.
    
    Returns
    -------
    tuple
        (is_valid : bool, score : float, explanation : str)
    """
    symbol = atom.GetSymbol()
    formal_charge = atom.GetFormalCharge()
    degree = atom.GetDegree()
    
    expected_max_valence = _get_expected_max_valence(symbol, formal_charge, strictness)
    
    if degree <= expected_max_valence:
        return True, 0.0, f"{symbol} valency correct ({degree}/{expected_max_valence})"
    else:
        excess = degree - expected_max_valence
        penalty = min(1.0, (excess ** 2) * 0.25)
        return False, penalty, f"{symbol} valency exceeded by {excess} (actual:{degree}, expected:{expected_max_valence})"


def check_aromaticity_detailed(atom, strictness='balanced'):
    """
    Check aromaticity consistency for atoms.
    
    Parameters
    ----------
    atom : rdkit.Chem.rdchem.Atom
        RDKit atom object.
    strictness : str, optional
        'conservative', 'balanced', or 'liberal'. Default 'balanced'.
    
    Returns
    -------
    tuple
        (is_valid : bool, score : float, explanation : str)
    """
    symbol = atom.GetSymbol()
    is_aromatic = atom.GetIsAromatic()
    is_in_ring = atom.IsInRing()
    
    # Allow aromatic atoms only if they're in rings
    if is_aromatic and not is_in_ring:
        return False, 0.3, f"{symbol} marked aromatic but not in ring"
    
    # For certain elements, check if aromatic form is reasonable
    aromatic_allowed = {'C', 'N', 'O', 'S'}
    if is_aromatic and symbol not in aromatic_allowed:
        return False, 0.2, f"{symbol} in aromatic context not typical"
    
    return True, 0.0, f"{symbol} aromaticity consistent"


def check_charge_detailed(atom, strictness='balanced'):
    """
    Check if atom's formal charge is chemically reasonable for its element.
    
    Parameters
    ----------
    atom : rdkit.Chem.rdchem.Atom
        RDKit atom object.
    strictness : str, optional
        'conservative', 'balanced', or 'liberal'. Default 'balanced'.
    
    Returns
    -------
    tuple
        (is_valid : bool, score : float, explanation : str)
    """
    symbol = atom.GetSymbol()
    formal_charge = atom.GetFormalCharge()
    
    # Define allowed charges for common elements
    allowed_charges = {
        'H': [-1, 0, 1],
        'C': [-1, 0, 1],      # Carbanion, neutral, carbocation/methyl
        'N': [-2, -1, 0, 1, 2],  # Various nitrogen oxidation states
        'O': [-1, 0, 1],      # Oxide, neutral, oxonium
        'F': [-1],            # Fluoride only
        'P': [-3, -2, -1, 0, 1, 2, 3, 4, 5],  # Phosphorus has many states
        'S': [-2, -1, 0, 1, 2, 4, 6],  # Sulfur oxidation states
        'Cl': [-1, 0, 1],     # Chlorine oxidation states
        'Br': [-1, 0, 1],     # Bromine oxidation states
        'I': [-1, 0, 1, 3, 5, 7]  # Iodine oxidation states
    }
    
    if symbol in allowed_charges:
        if formal_charge not in allowed_charges[symbol]:
            # Calculate distance from nearest allowed charge
            distances = [abs(formal_charge - allowed) for allowed in allowed_charges[symbol]]
            min_distance = min(distances)
            penalty = min(0.8, min_distance * 0.2)  # Cap at 0.8 for severe violations
            return False, penalty, f"{symbol} charge {formal_charge} not allowed (allowed: {allowed_charges[symbol]})"
        else:
            return True, 0.0, f"{symbol} charge {formal_charge} acceptable"
    else:
        # For elements not in our list, allow up to valence-based charges
        # Default to allowing charges that don't exceed typical valences
        if abs(formal_charge) <= 4:  # Reasonable limit for unknown elements
            return True, 0.0, f"{symbol} charge {formal_charge} within limits"
        else:
            penalty = min(0.8, (abs(formal_charge) - 4) * 0.2)
            return False, penalty, f"{symbol} charge {formal_charge} too extreme"


def _get_expected_max_valence(symbol, formal_charge, strictness='balanced'):
    """
    Helper function to get expected maximum valence based on element, formal charge, and strictness.
    
    Parameters
    ----------
    symbol : str
        Element symbol.
    formal_charge : int
        Formal charge of the atom.
    strictness : str
        'conservative', 'balanced', or 'liberal'.
    
    Returns
    -------
    int
        Expected maximum valence for the atom.
    """
    # Base valences for conservative strictness
    base_valences = {
        'H': 1, 'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1,
        'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2,
        'B': 3, 'Al': 3, 'Ga': 3, 'In': 3, 'Tl': 3,
        'C': 4, 'Si': 4, 'Ge': 4, 'Sn': 4, 'Pb': 4,
        'N': 3, 'P': 3, 'As': 3, 'Sb': 3, 'Bi': 3,
        'O': 2, 'S': 2, 'Se': 2, 'Te': 2, 'Po': 2,
        'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'At': 1
    }
    
    # Strictness-specific adjustments for sulfur and phosphorus
    if strictness == 'conservative':
        strictness_adjustments = {
            'S': 4,  # Conservative: sulfones not allowed
            'P': 4,  # Conservative: phosphates not allowed
        }
    elif strictness == 'balanced':
        strictness_adjustments = {
            'S': 6,  # Balanced: sulfones allowed
            'P': 5,  # Balanced: phosphates allowed
        }
    elif strictness == 'liberal':
        strictness_adjustments = {
            'S': 6,  # Liberal: sulfones allowed
            'P': 5,  # Liberal: phosphates allowed
        }
    else:
        strictness_adjustments = {}
    
    if symbol in base_valences:
        base_valence = base_valences[symbol]
        
        # Use strictness-specific adjustment for sulfur and phosphorus
        if symbol in strictness_adjustments:
            adjusted_valence = strictness_adjustments[symbol]
        else:
            # For other elements, adjust based on formal charge
            adjusted_valence = base_valence + formal_charge
            
        return max(1, adjusted_valence)  # Ensure at least 1
    else:
        # Default for unknown elements - allow up to 4 bonds typically
        return 4


def get_all_rules():
    """
    Return dictionary mapping rule names to their detailed check functions.
    
    Returns
    -------
    dict
        Dictionary of rule names to check functions.
    """
    return {
        'valency_check': check_valency_detailed,
        'aromaticity_check': check_aromaticity_detailed,
        'charge_check': check_charge_detailed,
        'ring_strain_check': lambda atom, strictness='balanced': (True, 0.0, "Handled at molecular level")  # Simplified placeholder
    }


if __name__ == "__main__":
    # Test the rules
    from rdkit import Chem
    
    # Test cases
    test_cases = [
        ("[NH4+]", "Nitrogen with +1 charge"),
        ("[C+]", "Carbon with +1 charge"),
        ("[C-]", "Carbon with -1 charge"),
        ("C", "Neutral carbon"),
        ("O", "Neutral oxygen"),
        ("[O-]", "Oxygen with -1 charge"),
        ("[OH3+]", "Oxygen with +1 charge"),
        ("N", "Neutral nitrogen"),
        ("[NH2-]", "Nitrogen with -1 charge"),
        ("[NH4+]", "Nitrogen with +1 charge"),
    ]
    
    for smiles, description in test_cases:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            atom = mol.GetAtomWithIdx(0)
            print(f"\n{description} ({smiles}):")
            print(f"  Valency: {check_valency_detailed(atom)}")
            print(f"  Charge: {check_charge_detailed(atom)}")
            print(f"  Aromaticity: {check_aromaticity_detailed(atom)}")