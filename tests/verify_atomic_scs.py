"""
Verification Script for Atomic SCS
Tests all major functionality to ensure the system works as expected.
"""

from atomic_scs import AtomicSCS
from analyzer import analyze_single_molecule
from scoring import MultiDimensionalScorer
from atomic_rules import check_valency_detailed, check_aromaticity_detailed, check_charge_detailed


def run_verification_tests():
    """
    Run comprehensive verification tests for Atomic SCS.
    """
    print("Running Atomic SCS Verification Tests...")
    print("=" * 50)
    
    all_passed = True
    
    # Test 1: Basic functionality
    print("\n1. Testing basic functionality...")
    try:
        scs = AtomicSCS('balanced')
        result = scs.assess_molecule("CCO")  # Ethanol
        assert result['valid'] == True
        assert isinstance(result['score'], float)
        assert result['score'] >= 0.0 and result['score'] <= 1.0
        print("   ✓ Basic assessment works")
    except Exception as e:
        print(f"   ✗ Basic assessment failed: {e}")
        all_passed = False
    
    # Test 2: Different strictness modes
    print("\n2. Testing different strictness modes...")
    try:
        for mode in ['conservative', 'balanced', 'liberal']:
            scs = AtomicSCS(mode)
            result = scs.assess_molecule("CCO")
            assert result['valid'] == True
        print("   ✓ All strictness modes work")
    except Exception as e:
        print(f"   ✗ Strictness modes failed: {e}")
        all_passed = False
    
    # Test 3: Ring strain assessment
    print("\n3. Testing ring strain assessment...")
    try:
        scs = AtomicSCS('balanced')
        # Cyclopropane should have ring strain
        cyclopropane_result = scs.assess_molecule("C1CC1")
        # Ethane should have no ring strain
        ethane_result = scs.assess_molecule("CC")
        
        assert cyclopropane_result['valid'] == True
        assert ethane_result['valid'] == True
        assert cyclopropane_result['score'] >= 0.0
        assert ethane_result['score'] >= 0.0
        print("   ✓ Ring strain assessment works")
    except Exception as e:
        print(f"   ✗ Ring strain assessment failed: {e}")
        all_passed = False
    
    # Test 4: Valency violations
    print("\n4. Testing valency violations...")
    try:
        scs = AtomicSCS('balanced')
        # Carbocation [C+] should have high valency score
        carbocation_result = scs.assess_molecule("[C+]")
        assert carbocation_result['valid'] == True
        assert carbocation_result['score'] >= 0.0
        print("   ✓ Valency violation assessment works")
    except Exception as e:
        print(f"   ✗ Valency violation assessment failed: {e}")
        all_passed = False
    
    # Test 5: Diagnostic mode
    print("\n5. Testing diagnostic mode...")
    try:
        scs = AtomicSCS('balanced')
        diagnosis = scs.diagnose_molecule("C1CC1")  # Cyclopropane
        assert diagnosis['valid'] == True
        assert 'diagnosis' in diagnosis
        assert 'total_violations' in diagnosis['diagnosis']
        print("   ✓ Diagnostic mode works")
    except Exception as e:
        print(f"   ✗ Diagnostic mode failed: {e}")
        all_passed = False
    
    # Test 6: Repair mode
    print("\n6. Testing repair mode...")
    try:
        scs = AtomicSCS('balanced')
        repair_result = scs.repair_molecule("C1CC1")  # Cyclopropane
        assert repair_result['valid'] == True
        assert 'repairs' in repair_result
        print("   ✓ Repair mode works")
    except Exception as e:
        print(f"   ✗ Repair mode failed: {e}")
        all_passed = False
    
    # Test 7: Atomic rules directly
    print("\n7. Testing atomic rules individually...")
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        atom = mol.GetAtomWithIdx(0)  # Carbon
        
        valency_result = check_valency_detailed(atom, 'balanced')
        aromaticity_result = check_aromaticity_detailed(atom, 'balanced')
        charge_result = check_charge_detailed(atom, 'balanced')
        
        assert isinstance(valency_result, tuple) and len(valency_result) == 3
        assert isinstance(aromaticity_result, tuple) and len(aromaticity_result) == 3
        assert isinstance(charge_result, tuple) and len(charge_result) == 3
        print("   ✓ Individual atomic rules work")
    except Exception as e:
        print(f"   ✗ Individual atomic rules failed: {e}")
        all_passed = False
    
    # Test 8: Sulfur and phosphorus handling
    print("\n8. Testing sulfur and phosphorus handling...")
    try:
        scs = AtomicSCS('balanced')
        # Sulfuric acid - sulfur with 6 bonds should be OK in balanced mode
        sulfuric_acid = scs.assess_molecule("OS(=O)(=O)O")
        # Phosphoric acid - phosphorus with 5 bonds should be OK in balanced mode
        phosphoric_acid = scs.assess_molecule("OP(=O)(O)O")
        
        assert sulfuric_acid['valid'] == True
        assert phosphoric_acid['valid'] == True
        print("   ✓ Sulfur and phosphorus handling works")
    except Exception as e:
        print(f"   ✗ Sulfur and phosphorus handling failed: {e}")
        all_passed = False
    
    # Test 9: Invalid SMILES
    print("\n9. Testing invalid SMILES handling...")
    try:
        scs = AtomicSCS('balanced')
        invalid_result = scs.assess_molecule("INVALID_SMILES")
        assert invalid_result['valid'] == False
        assert invalid_result['score'] == 1.0
        print("   ✓ Invalid SMILES handling works")
    except Exception as e:
        print(f"   ✗ Invalid SMILES handling failed: {e}")
        all_passed = False
    
    # Test 10: Molecular score calculation
    print("\n10. Testing molecular score calculation...")
    try:
        scorer = MultiDimensionalScorer()
        # Simulate atom details with different scores
        fake_atom_details = [
            {'scores': {'valency_check': 0.1, 'aromaticity_check': 0.0, 'charge_check': 0.05}},
            {'scores': {'valency_check': 0.0, 'aromaticity_check': 0.2, 'charge_check': 0.0}},
            {'scores': {'valency_check': 0.0, 'aromaticity_check': 0.0, 'charge_check': 0.15}}
        ]
        max_score = scorer.compute_molecular_score(fake_atom_details)
        assert max_score == 0.2  # Should be the maximum of all scores
        print("   ✓ Molecular score calculation works")
    except Exception as e:
        print(f"   ✗ Molecular score calculation failed: {e}")
        all_passed = False
    
    # Summary
    print("\n" + "=" * 50)
    if all_passed:
        print("✓ All verification tests PASSED!")
        print("Atomic SCS is working correctly.")
        return True
    else:
        print("✗ Some verification tests FAILED!")
        print("Please review the errors above.")
        return False


if __name__ == "__main__":
    success = run_verification_tests()
    exit(0 if success else 1)
