#!/usr/bin/env python3
"""
Benchmark script for Atomic-SCS.
Usage: python benchmark.py <smiles_file>
Example: python benchmark.py normal_100.smi
"""

import sys
import time
from atomic_scs import AtomicSCS

def benchmark(smiles_file, strictness='balanced'):
    """Run performance benchmark on a SMILES file."""
    scs = AtomicSCS(strictness=strictness)
    
    # Read SMILES
    with open(smiles_file, 'r') as f:
        smiles_list = [line.strip().split()[0] for line in f if line.strip()]
    
    if not smiles_list:
        print("Error: No SMILES found.")
        return
    
    print(f"Benchmarking {len(smiles_list)} molecules...")
    
    # Warm-up (optional, to avoid first-call overhead)
    scs.assess_molecule(smiles_list[0])
    
    # Measure time
    start = time.perf_counter()
    for smiles in smiles_list:
        scs.assess_molecule(smiles)
    end = time.perf_counter()
    
    total = end - start
    avg = total / len(smiles_list)
    
    print(f"Total time: {total:.4f} seconds")
    print(f"Average per molecule: {avg:.6f} seconds")
    print(f"Molecules per second: {1/avg:.1f}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python benchmark.py <smiles_file>")
        sys.exit(1)
    benchmark(sys.argv[1])