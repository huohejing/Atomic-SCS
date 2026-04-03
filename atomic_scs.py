"""
Atomic SCS (Structure Compliance System)
Main module for atomic-level molecular compliance assessment
"""

from rdkit import Chem
from analyzer import MoleculeAnalyzer, analyze_single_molecule
from scoring import MultiDimensionalScorer
import argparse
import sys
import os
import csv
import json


class AtomicSCS:
    def __init__(self, strictness='balanced'):
        """
        Initialize AtomicSCS with specified strictness mode
        
        Args:
            strictness: 'conservative', 'balanced', or 'liberal'
        """
        self.strictness = strictness.lower()
        if self.strictness not in ['conservative', 'balanced', 'liberal']:
            raise ValueError("Strictness must be 'conservative', 'balanced', or 'liberal'")
        
        self.analyzer = MoleculeAnalyzer()
        self.scorer = MultiDimensionalScorer()
    
    def assess_molecule(self, smiles):
        """
        Assess a single molecule for compliance
        
        Args:
            smiles: SMILES string of the molecule
        
        Returns:
            dict: Assessment results
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                'valid': False,
                'smiles': smiles,
                'error': 'Invalid SMILES',
                'score': 1.0,
                'compliant': False,
                'confidence_level': 'very_poor'
            }
        
        # Analyze the molecule
        analysis_result = self.analyzer.analyze_molecule(mol, self.strictness)
        
        if not analysis_result['valid']:
            return {
                'valid': False,
                'smiles': smiles,
                'error': analysis_result.get('error', 'Analysis failed'),
                'score': 1.0,
                'compliant': False,
                'confidence_level': 'very_poor'
            }
        
        # Calculate molecular score as the maximum of all atom scores across all rules
        molecular_score = self.scorer.compute_molecular_score(analysis_result['atom_details'])
        
        # Assess compliance based on strictness-specific threshold
        compliant = self.scorer.assess_compliance(molecular_score, self.strictness)
        confidence_level = self.scorer.get_confidence_level(molecular_score)
        
        return {
            'valid': True,
            'smiles': smiles,
            'score': molecular_score,
            'compliant': compliant,
            'confidence_level': confidence_level,
            'atom_details': analysis_result['atom_details'],
            'num_atoms': analysis_result['num_atoms']
        }
    
    def diagnose_molecule(self, smiles):
        """
        Provide detailed diagnosis of compliance issues
        """
        assessment = self.assess_molecule(smiles)
        
        if not assessment['valid']:
            return assessment
        
        # Generate repair guidance
        guidance = self.scorer.generate_repair_guidance(assessment)
        
        assessment['diagnosis'] = {
            'total_violations': len(guidance),
            'violations_by_severity': self._count_by_severity(guidance),
            'repair_guidance': guidance
        }
        
        return assessment
    
    def repair_molecule(self, smiles):
        """
        Suggest repairs for compliance issues
        """
        diagnosis = self.diagnose_molecule(smiles)
        
        if not diagnosis['valid']:
            return diagnosis
        
        # Return prioritized repair suggestions
        diagnosis['repairs'] = diagnosis['diagnosis']['repair_guidance']
        
        return diagnosis
    
    def _count_by_severity(self, guidance):
        """
        Count violations by severity level
        """
        counts = {}
        for item in guidance:
            severity = item['severity']
            counts[severity] = counts.get(severity, 0) + 1
        return counts


def write_text_output(out_f, mode, verbose):
    """Write text (TSV) format headers"""
    if mode == 'assess':
        out_f.write("SMILES\tSCORE\tCOMPLIANT\tCONFIDENCE_LEVEL\n")
    elif mode == 'diagnose':
        if verbose:
            out_f.write("SMILES\tATOM_INDEX\tELEMENT\tRULE_TYPE\tSCORE\tEXPLANATION\n")
        else:
            out_f.write("SMILES\tSCORE\tCOMPLIANT\tVIOLATIONS_COUNT\tSEVERITY_BREAKDOWN\n")
    elif mode == 'repair':
        if verbose:
            out_f.write("SMILES\tATOM_INDEX\tELEMENT\tSEVERITY\tSCORE\tSUGGESTION\n")
        else:
            out_f.write("SMILES\tSCORE\tCOMPLIANT\tREPAIR_SUGGESTIONS_COUNT\tTOP_REPAIRS\n")


def write_csv_output(csv_writer, mode, verbose):
    """Write CSV format headers"""
    if mode == 'assess':
        csv_writer.writerow(["SMILES", "SCORE", "COMPLIANT", "CONFIDENCE_LEVEL"])
    elif mode == 'diagnose':
        if verbose:
            csv_writer.writerow(["SMILES", "ATOM_INDEX", "ELEMENT", "RULE_TYPE", "SCORE", "EXPLANATION"])
        else:
            csv_writer.writerow(["SMILES", "SCORE", "COMPLIANT", "VIOLATIONS_COUNT", "SEVERITY_BREAKDOWN"])
    elif mode == 'repair':
        if verbose:
            csv_writer.writerow(["SMILES", "ATOM_INDEX", "ELEMENT", "SEVERITY", "SCORE", "SUGGESTION"])
        else:
            csv_writer.writerow(["SMILES", "SCORE", "COMPLIANT", "REPAIR_SUGGESTIONS_COUNT", "TOP_REPAIRS"])


def process_file(input_file, output_file=None, strictness='balanced', mode='assess', verbose=False, output_format='text'):
    """
    Process a file containing SMILES strings
    
    Args:
        input_file: Path to input file with SMILES
        output_file: Path to output file (None for stdout)
        strictness: Strictness level
        mode: Operation mode
        verbose: Whether to output detailed atom-level information
        output_format: Output format ('text', 'csv', or 'json')
    """
    scs = AtomicSCS(strictness)
    
    # Determine output destination
    if output_file:
        out_f = open(output_file, 'w', newline='')  # newline='' for csv compatibility
    else:
        out_f = sys.stdout

    try:
        results_list = []  # For storing results when using JSON format
        
        # Handle JSON format separately since it needs all results at once
        if output_format == 'json':
            # Process all molecules first
            with open(input_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    smiles = line.strip().split()[0] if line.strip() else ""  # Take first column if tab-separated
                    
                    if not smiles:
                        continue
                    
                    try:
                        if mode == 'assess':
                            result = scs.assess_molecule(smiles)
                            json_result = {
                                'smiles': smiles,
                                'score': round(result['score'], 6),
                                'compliant': result['compliant'],
                                'confidence_level': result['confidence_level']
                            }
                            results_list.append(json_result)
                        
                        elif mode == 'diagnose':
                            result = scs.diagnose_molecule(smiles)
                            if verbose and result['valid']:
                                # Detailed per-atom diagnosis
                                has_violations = False
                                for atom_data in result['atom_details']:
                                    if atom_data['violations']:
                                        has_violations = True
                                        for violation in atom_data['violations']:
                                            json_result = {
                                                'smiles': smiles,
                                                'atom_index': atom_data['index'],
                                                'element': atom_data['element'],
                                                'rule_type': violation['rule'],
                                                'score': round(violation['score'], 6),
                                                'explanation': violation['explanation']
                                            }
                                            results_list.append(json_result)
                                
                                # If no violations found, add an OK status record
                                if not has_violations:
                                    json_result = {
                                        'smiles': smiles,
                                        'status': 'OK',
                                        'message': 'No violations detected'
                                    }
                                    results_list.append(json_result)
                            else:
                                # Summary diagnosis
                                violations_count = result.get('diagnosis', {}).get('total_violations', 0)
                                severity_breakdown = result.get('diagnosis', {}).get('violations_by_severity', {})
                                json_result = {
                                    'smiles': smiles,
                                    'score': round(result['score'], 6),
                                    'compliant': result['compliant'],
                                    'violations_count': violations_count,
                                    'severity_breakdown': severity_breakdown
                                }
                                results_list.append(json_result)
                        
                        elif mode == 'repair':
                            result = scs.repair_molecule(smiles)
                            if verbose and result['valid']:
                                # Detailed repair suggestions
                                has_repairs = False
                                for repair in result.get('repairs', []):
                                    has_repairs = True
                                    json_result = {
                                        'smiles': smiles,
                                        'atom_index': repair['atom_index'],
                                        'element': repair['element'],
                                        'severity': repair['severity'],
                                        'score': round(repair['score'], 6),
                                        'suggestion': repair['suggestion']
                                    }
                                    results_list.append(json_result)
                                
                                # If no repairs needed, add an OK status record
                                if not has_repairs:
                                    json_result = {
                                        'smiles': smiles,
                                        'status': 'OK',
                                        'message': 'No repairs needed'
                                    }
                                    results_list.append(json_result)
                            else:
                                # Summary repair info
                                repairs_count = len(result.get('repairs', []))
                                top_repairs = [r['suggestion'] for r in result.get('repairs', [])[:3]]  # Top 3 repairs
                                json_result = {
                                    'smiles': smiles,
                                    'score': round(result['score'], 6),
                                    'compliant': result['compliant'],
                                    'repair_suggestions_count': repairs_count,
                                    'top_repairs': top_repairs
                                }
                                results_list.append(json_result)
                    
                    except Exception as e:
                        error_msg = f"Error processing line {line_num}: {str(e)}"
                        print(error_msg, file=sys.stderr)
                        # Add error entry to results
                        if mode == 'assess':
                            results_list.append({
                                'smiles': smiles,
                                'score': 1.0,
                                'compliant': False,
                                'confidence_level': 'ERROR',
                                'error': str(e)
                            })
                        elif mode == 'diagnose':
                            if verbose:
                                results_list.append({
                                    'smiles': smiles,
                                    'atom_index': -1,
                                    'element': 'ERROR',
                                    'rule_type': 'ERROR',
                                    'score': 1.0,
                                    'explanation': str(e)
                                })
                            else:
                                results_list.append({
                                    'smiles': smiles,
                                    'score': 1.0,
                                    'compliant': False,
                                    'violations_count': 0,
                                    'severity_breakdown': {},
                                    'error': str(e)
                                })
                        elif mode == 'repair':
                            if verbose:
                                results_list.append({
                                    'smiles': smiles,
                                    'atom_index': -1,
                                    'element': 'ERROR',
                                    'severity': 'ERROR',
                                    'score': 1.0,
                                    'suggestion': str(e)
                                })
                            else:
                                results_list.append({
                                    'smiles': smiles,
                                    'score': 1.0,
                                    'compliant': False,
                                    'repair_suggestions_count': 0,
                                    'top_repairs': [],
                                    'error': str(e)
                                })
            
            # Write all results as JSON array
            json.dump(results_list, out_f, indent=2)
        
        else:  # Handle text and csv formats
            if output_format == 'csv':
                csv_writer = csv.writer(out_f)
                write_csv_output(csv_writer, mode, verbose)
            else:  # text format (default TSV)
                write_text_output(out_f, mode, verbose)
            
            with open(input_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    smiles = line.strip().split()[0] if line.strip() else ""  # Take first column if tab-separated
                    
                    if not smiles:
                        continue
                    
                    try:
                        if mode == 'assess':
                            result = scs.assess_molecule(smiles)
                            
                            if output_format == 'csv':
                                csv_writer.writerow([
                                    smiles,
                                    f"{result['score']:.6f}",
                                    result['compliant'],
                                    result['confidence_level']
                                ])
                            else:  # text (TSV)
                                out_f.write(f"{smiles}\t{result['score']:.6f}\t{result['compliant']}\t{result['confidence_level']}\n")
                        
                        elif mode == 'diagnose':
                            if verbose:
                                # Output detailed atom-level information
                                result = scs.diagnose_molecule(smiles)
                                if result['valid']:
                                    has_violations = False
                                    for atom_data in result['atom_details']:
                                        if atom_data['violations']:
                                            has_violations = True
                                            for violation in atom_data['violations']:
                                                if output_format == 'csv':
                                                    csv_writer.writerow([
                                                        smiles,
                                                        atom_data['index'],
                                                        atom_data['element'],
                                                        violation['rule'],
                                                        f"{violation['score']:.6f}",
                                                        violation['explanation']
                                                    ])
                                                else:  # text (TSV)
                                                    out_f.write(f"{smiles}\t{atom_data['index']}\t{atom_data['element']}\t"
                                                               f"{violation['rule']}\t{violation['score']:.6f}\t{violation['explanation']}\n")
                                    
                                    # If no violations were found, still output a message
                                    if not has_violations:
                                        if output_format == 'csv':
                                            csv_writer.writerow([smiles, "OK", "No violations detected"])
                                        else:  # text (TSV)
                                            out_f.write(f"{smiles}\tOK\tNo violations detected\n")
                            else:
                                # Output summary information
                                result = scs.diagnose_molecule(smiles)
                                violations_count = result.get('diagnosis', {}).get('total_violations', 0)
                                severity_breakdown = str(result.get('diagnosis', {}).get('violations_by_severity', {}))
                                
                                if output_format == 'csv':
                                    csv_writer.writerow([
                                        smiles,
                                        f"{result['score']:.6f}",
                                        result['compliant'],
                                        violations_count,
                                        severity_breakdown
                                    ])
                                else:  # text (TSV)
                                    out_f.write(f"{smiles}\t{result['score']:.6f}\t{result['compliant']}\t{violations_count}\t{severity_breakdown}\n")
                        
                        elif mode == 'repair':
                            if verbose:
                                # Output detailed repair suggestions
                                result = scs.repair_molecule(smiles)
                                if result['valid']:
                                    has_repairs = False
                                    for repair in result.get('repairs', []):
                                        has_repairs = True
                                        if output_format == 'csv':
                                            csv_writer.writerow([
                                                smiles,
                                                repair['atom_index'],
                                                repair['element'],
                                                repair['severity'],
                                                f"{repair['score']:.6f}",
                                                repair['suggestion']
                                            ])
                                        else:  # text (TSV)
                                            out_f.write(f"{smiles}\t{repair['atom_index']}\t{repair['element']}\t"
                                                       f"{repair['severity']}\t{repair['score']:.6f}\t{repair['suggestion']}\n")
                                    
                                    # If no repairs were needed, still output a message
                                    if not has_repairs:
                                        if output_format == 'csv':
                                            csv_writer.writerow([smiles, "OK", "No repairs needed"])
                                        else:  # text (TSV)
                                            out_f.write(f"{smiles}\tOK\tNo repairs needed\n")
                            else:
                                # Output summary information
                                result = scs.repair_molecule(smiles)
                                repairs_count = len(result.get('repairs', []))
                                top_repairs = [r['suggestion'] for r in result.get('repairs', [])[:3]]  # Top 3 repairs
                                top_repairs_str = "; ".join(top_repairs)
                                
                                if output_format == 'csv':
                                    csv_writer.writerow([
                                        smiles,
                                        f"{result['score']:.6f}",
                                        result['compliant'],
                                        repairs_count,
                                        top_repairs_str
                                    ])
                                else:  # text (TSV)
                                    out_f.write(f"{smiles}\t{result['score']:.6f}\t{result['compliant']}\t{repairs_count}\t{top_repairs_str}\n")
                    
                    except Exception as e:
                        error_msg = f"Error processing line {line_num}: {str(e)}"
                        if output_format == 'csv':
                            if mode == 'assess':
                                csv_writer.writerow([smiles, "1.0", "False", "ERROR"])
                            elif mode == 'diagnose':
                                if verbose:
                                    csv_writer.writerow([smiles, "-1", "ERROR", "ERROR", "1.0", str(e)])
                                else:
                                    csv_writer.writerow([smiles, "1.0", "False", "0", "{}"])
                            elif mode == 'repair':
                                if verbose:
                                    csv_writer.writerow([smiles, "-1", "ERROR", "ERROR", "1.0", str(e)])
                                else:
                                    csv_writer.writerow([smiles, "1.0", "False", "0", ""])
                        else:  # text (TSV)
                            if mode == 'assess':
                                out_f.write(f"{smiles}\t1.0\tFalse\tERROR\n")
                            elif mode == 'diagnose':
                                if verbose:
                                    out_f.write(f"{smiles}\t-1\tERROR\tERROR\t1.0\t{str(e)}\n")
                                else:
                                    out_f.write(f"{smiles}\t1.0\tFalse\t0\t{{}}\n")
                            elif mode == 'repair':
                                if verbose:
                                    out_f.write(f"{smiles}\t-1\tERROR\tERROR\t1.0\t{str(e)}\n")
                                else:
                                    out_f.write(f"{smiles}\t1.0\tFalse\t0\t\n")
                        print(error_msg, file=sys.stderr)
    
    finally:
        if output_file:
            out_f.close()


def main():
    parser = argparse.ArgumentParser(description="Atomic Structure Compliance System")
    parser.add_argument("input_file", help="Input file with SMILES strings")
    parser.add_argument("--output", "-o", help="Output file (default: stdout)")
    parser.add_argument("--strictness", "-s", choices=['conservative', 'balanced', 'liberal'], 
                       default='balanced', help="Strictness level (default: balanced)")
    parser.add_argument("--mode", "-m", choices=['assess', 'diagnose', 'repair'], 
                       default='assess', help="Operation mode (default: assess)")
    parser.add_argument("--verbose", "-v", action='store_true', 
                       help="Output detailed atom-level information")
    parser.add_argument("--format", "-f", choices=['text', 'csv', 'json'], default='text',
                        help="Output format (default: text)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    process_file(args.input_file, args.output, args.strictness, args.mode, args.verbose, args.format)


if __name__ == "__main__":
    # Example usage
    if len(sys.argv) == 1:
        # Run examples when called without arguments
        scs = AtomicSCS('balanced')
        
        examples = [
            "CCO",           # Ethanol
            "C1CC1",         # Cyclopropane (should have ring strain)
            "[NH4+]",        # Ammonium ion
            "C[CH2+]",       # Carbocation
            "CC(=O)O",       # Acetic acid
            "OS(=O)(=O)O",   # Sulfuric acid (S=6)
            "OP(=O)(O)O",    # Phosphoric acid (P=5)
        ]
        
        print("Example assessments:")
        for smiles in examples:
            result = scs.assess_molecule(smiles)
            print(f"{smiles}: score={result['score']:.3f}, compliant={result['compliant']}, "
                  f"confidence={result['confidence_level']}")
        
        print("\nDetailed diagnosis for cyclopropane:")
        diagnosis = scs.diagnose_molecule("C1CC1")
        if diagnosis['valid']:
            print(f"Score: {diagnosis['score']:.3f}")
            print(f"Violations: {diagnosis['diagnosis']['total_violations']}")
            for violation in diagnosis['diagnosis']['repair_guidance'][:3]:
                print(f"  {violation['element']}@{violation['atom_index']}: {violation['explanation']}")
        
        print("\nTest for sulfuric acid (should handle S=6 in balanced mode):")
        sulfuric_acid = scs.assess_molecule("OS(=O)(=O)O")
        print(f"Sulfuric acid: score={sulfuric_acid['score']:.3f}, compliant={sulfuric_acid['compliant']}")
    else:
        # Run command-line interface
        main()