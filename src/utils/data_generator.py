"""
Universal Data Generator for Genomic Pattern Matching Analysis.

Handles:
1. Synthetic DNA Generation (10^4 to 10^7 bp)
2. Real Data Loading (FASTA/FASTQ parsing)
3. Pattern Extraction (Short, Medium, Long)
4. Mutation Injection (for Approximate Matching testing)

Author: Team Subproblem Solvers
"""

import random
import os

def load_genome_from_file(filepath):
    """
    Loads a DNA sequence from a FASTA or raw text file.
    
    Args:
        filepath (str): Path to the .fasta or .txt file.
        
    Returns:
        str: The complete DNA sequence as a single string (uppercase).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    sequence = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip FASTA headers (lines starting with >)
            if line.startswith(">") or line.startswith("@"):
                continue
            sequence.append(line.upper())
    
    # Join all lines into one massive string
    full_sequence = "".join(sequence)
    
    # rudimentary cleaning to ensure only DNA chars exist
    # (Optional: remove N or other artifacts if necessary)
    return full_sequence


def generate_synthetic_dna(length, seed=42):
    """
    Generates a random DNA sequence of exact length.
    
    Args:
        length (int): Number of base pairs (e.g., 10**7).
        seed (int): Random seed for reproducibility.
        
    Returns:
        str: Random DNA sequence.
    """
    random.seed(seed)
    return ''.join(random.choices(['A', 'T', 'C', 'G'], k=length))


def generate_patterns_from_text(text, length, count=1, seed=None):
    """
    Extracts real patterns that actually exist in the text.
    Use this for Exact Matching Benchmarks (KMP, Boyer-Moore).
    
    Args:
        text (str): The large DNA sequence.
        length (int): Length of the motif (e.g., 10, 20, 50).
        count (int): How many different patterns to extract.
        
    Returns:
        list[str]: A list of patterns guaranteed to be in the text.
    """
    if seed:
        random.seed(seed)
        
    patterns = []
    text_len = len(text)
    
    if length > text_len:
        raise ValueError("Pattern length cannot be larger than text length.")

    for _ in range(count):
        # Pick a random starting point
        start_idx = random.randint(0, text_len - length)
        # Extract the slice
        p = text[start_idx : start_idx + length]
        patterns.append(p)
        
    return patterns

def generate_mutated_pattern(original_pattern, k, mutation_type="mix", seed=None):
    """
    Takes a valid pattern and introduces 'k' errors to simulate real sequencing data.
    
    Args:
        original_pattern (str): The clean DNA motif.
        k (int): Number of edit operations (errors) to introduce.
        mutation_type (str): 
            - "sub": Only substitutions (Good for Shift-Or basic).
            - "mix": Mix of insertions, deletions, substitutions (Good for Levenshtein).
        seed (int): For reproducibility.
        
    Returns:
        str: The mutated pattern.
    """
    if seed:
        random.seed(seed)
        
    pattern_list = list(original_pattern)
    bases = ['A', 'T', 'C', 'G']
    
    for _ in range(k):
        # Decide what kind of mutation to do
        if mutation_type == "sub":
            op = "sub"
        else:
            # Randomly choose between substitution, insertion, or deletion
            op = random.choice(["sub", "ins", "del"])
        
        # 1. SUBSTITUTION (Change a char)
        if op == "sub" and len(pattern_list) > 0:
            idx = random.randint(0, len(pattern_list) - 1)
            current_base = pattern_list[idx]
            choices = [b for b in bases if b != current_base]
            pattern_list[idx] = random.choice(choices)
            
        # 2. INSERTION (Add a char)
        elif op == "ins":
            idx = random.randint(0, len(pattern_list)) # Can insert at end
            pattern_list.insert(idx, random.choice(bases))
            
        # 3. DELETION (Remove a char)
        elif op == "del" and len(pattern_list) > 0:
            idx = random.randint(0, len(pattern_list) - 1)
            pattern_list.pop(idx)
            
    return "".join(pattern_list)

# --- QUICK TEST BLOCK (Runs only if you execute this file directly) ---
if __name__ == "__main__":
    print("--- Testing Universal Data Generator ---")
    
    # 1. Test Scale (10^6)
    print("Generating 1 Million Base Pairs...")
    large_dna = generate_synthetic_dna(10**6)
    print(f"Length: {len(large_dna)}")
    print(f"Preview: {large_dna[:50]}...")
    
    # 2. Test Pattern Extraction (Medium Length)
    print("\nExtracting Pattern (Length 20)...")
    patterns = generate_patterns_from_text(large_dna, 20, count=1)
    print(f"Pattern: {patterns[0]}")
    
    # 3. Test Mutation (k=3)
    print("\nCreating Mutated Version (k=3)...")
    mutated = generate_mutated_pattern(patterns[0], 3)
    print(f"Original: {patterns[0]}")
    print(f"Mutated : {mutated}")
    
    # Verify differences
    diff_count = sum(1 for a, b in zip(patterns[0], mutated) if a != b)
    print(f"Actual Differences: {diff_count}")