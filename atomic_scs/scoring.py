"""
Scoring Module
Implements multi-dimensional scoring with configurable weights
"""

class MultiDimensionalScorer:
    def __init__(self):
        # Define severity levels for repair suggestions
        self.severity_map = {
            'critical': 5,
            'high': 4,
            'medium': 3,
            'low': 2,
            'info': 1
        }
    
    def compute_molecular_score(self, atom_details):
        """
        Calculate molecular score as the maximum of all atom scores across all rules.
        This ensures severe violations are not diluted by compliant atoms.
        
        Parameters
        ----------
        atom_details : list
            List of dictionaries containing atom-level analysis results.
        
        Returns
        -------
        float
            Maximum score across all atoms and all rules.
        """
        max_score = 0.0
        for atom in atom_details:
            for rule_score in atom['scores'].values():
                max_score = max(max_score, rule_score)
        return max_score
    
    def get_compliance_threshold(self, strictness):
        """
        Get compliance threshold based on strictness level.
        
        Parameters
        ----------
        strictness : str
            'conservative', 'balanced', or 'liberal'.
        
        Returns
        -------
        float
            Threshold value for compliance assessment.
        """
        thresholds = {
            'conservative': 0.1,   # Lower threshold for conservative mode (more strict)
            'balanced': 0.3,       # Standard threshold
            'liberal': 0.5         # Higher threshold for liberal mode (less strict)
        }
        return thresholds.get(strictness, 0.3)
    
    def assess_compliance(self, normalized_score, strictness):
        """
        Assess whether the molecule complies with standards based on normalized score.
        
        Parameters
        ----------
        normalized_score : float
            The computed molecular score.
        strictness : str
            'conservative', 'balanced', or 'liberal'.
        
        Returns
        -------
        bool
            True if the molecule is compliant, False otherwise.
        """
        threshold = self.get_compliance_threshold(strictness)
        return normalized_score <= threshold
    
    def get_confidence_level(self, normalized_score):
        """
        Get confidence level based on score.
        
        Parameters
        ----------
        normalized_score : float
            The computed molecular score.
        
        Returns
        -------
        str
            Confidence level ('very_good', 'good', 'fair', 'poor', 'very_poor').
        """
        if normalized_score <= 0.05:
            return 'very_good'
        elif normalized_score <= 0.15:
            return 'good'
        elif normalized_score <= 0.3:
            return 'fair'
        elif normalized_score <= 0.5:
            return 'poor'
        else:
            return 'very_poor'
    
    def generate_repair_guidance(self, analysis_result):
        """
        Generate repair guidance based on violations found.
        Uses severity and raw score for sorting (no additional weighting).
        
        Parameters
        ----------
        analysis_result : dict
            Result from molecule analysis containing atom details and violations.
        
        Returns
        -------
        list
            Sorted list of repair suggestions ordered by severity and score.
        """
        if not analysis_result.get('valid', False):
            return [{'issue': 'Invalid molecule', 'suggestion': 'Check SMILES syntax'}]
        
        violations = []
        
        for atom_data in analysis_result.get('atom_details', []):
            for violation in atom_data.get('violations', []):
                severity = self._classify_severity(violation['rule'])
                
                violations.append({
                    'atom_index': atom_data['index'],
                    'element': atom_data['element'],
                    'rule_type': violation['rule'],
                    'severity': severity,
                    'score': violation['score'],
                    'explanation': violation['explanation'],
                    'suggestion': self._generate_suggestion(violation['rule'], atom_data)
                })
        
        # Sort by severity first, then by raw score (higher score = more severe issue)
        sorted_violations = sorted(
            violations,
            key=lambda x: (
                -self.severity_map.get(x['severity'], 0),
                -x['score']  # Sort by raw score (higher first)
            )
        )
        
        return sorted_violations
    
    def _classify_severity(self, rule_name):
        """
        Classify severity level based on rule type.
        
        Parameters
        ----------
        rule_name : str
            Name of the rule that was violated.
        
        Returns
        -------
        str
            Severity level ('critical', 'high', 'medium', 'low').
        """
        critical_rules = ['valency_check']
        high_rules = ['charge_check']
        medium_rules = ['aromaticity_check']
        
        if rule_name in critical_rules:
            return 'critical'
        elif rule_name in high_rules:
            return 'high'
        elif rule_name in medium_rules:
            return 'medium'
        else:
            return 'low'
    
    def _generate_suggestion(self, rule_name, atom_data):
        """
        Generate repair suggestion based on rule violation.
        
        Parameters
        ----------
        rule_name : str
            Name of the violated rule.
        atom_data : dict
            Information about the atom that violated the rule.
        
        Returns
        -------
        str
            Suggestion for how to fix the violation.
        """
        suggestions = {
            'valency_check': f"Reduce number of bonds to {atom_data['element']} atom. "
                            f"Current degree: {atom_data['degree']}.",
            'charge_check': f"Adjust formal charge on {atom_data['element']} atom. "
                           f"Current charge: {atom_data['formal_charge']}.",
            'aromaticity_check': f"Review aromaticity assignment for {atom_data['element']} atom.",
            'ring_strain_check': f"Consider larger ring size to reduce strain."
        }
        
        return suggestions.get(rule_name, f"Review {rule_name} for atom {atom_data['element']}")


def calculate_molecular_score(atom_details, method='max'):
    """
    Calculate overall molecular score from atom-level details.
    Deprecated: Use MultiDimensionalScorer.compute_molecular_score instead.
    
    Parameters
    ----------
    atom_details : list
        List of atom-level analysis results.
    method : str, optional
        Aggregation method ('max' only supported). Default 'max'.
    
    Returns
    -------
    float
        Molecular score.
    """
    scorer = MultiDimensionalScorer()
    return scorer.compute_molecular_score(atom_details)


if __name__ == "__main__":
    # Test the scoring functionality
    from analyzer import analyze_single_molecule
    
    test_molecules = [
        "C",           # Normal
        "C1CC1",       # Small ring
        "[NH4+]",      # Charged
        "C[CH2+]",     # Carbocation
        "C1CCCCC1",    # Normal ring
        "OS(=O)(=O)O", # Sulfuric acid (S=6)
        "OP(=O)(O)O",  # Phosphoric acid (P=5)
    ]
    
    for smiles in test_molecules:
        print(f"\nTesting: {smiles}")
        analysis = analyze_single_molecule(smiles)
        if analysis['valid']:
            score = calculate_molecular_score(analysis['atom_details'])
            print(f"  Molecular score (max): {score:.3f}")
            
            scorer = MultiDimensionalScorer()
            guidance = scorer.generate_repair_guidance(analysis)
            if guidance:
                print(f"  Violations found: {len(guidance)}")
                for g in guidance[:2]:  # Show first 2
                    print(f"    {g['element']}@{g['atom_index']}: {g['explanation']}")
            else:
                print("  No violations found")
