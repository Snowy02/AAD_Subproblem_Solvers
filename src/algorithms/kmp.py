"""
Knuth-Morris-Pratt (KMP) Algorithm Implementation.

This module provides an efficient exact pattern matching algorithm
that preprocesses the pattern to construct a Longest Prefix-Suffix (LPS) 
table, allowing it to skip unnecessary comparisons in the text.
"""

def compute_lps_table(pattern: str) -> list[int]:
    """
    Computes the LPS (Longest Proper Prefix which is also a Suffix) table.

    Args:
        pattern (str): The search pattern (DNA motif).

    Returns:
        list[int]: The LPS array where lps[i] is the length of the longest
                   proper prefix of pattern[0...i] that is also a suffix.
    """
    length = len(pattern)
    lps = [0] * length
    j = 0  # Length of the previous longest prefix suffix
    i = 1

    while i < length:
        if pattern[i] == pattern[j]:
            j += 1
            lps[i] = j
            i += 1
        else:
            if j != 0:
                # Backtrack to the previous longest prefix-suffix.
                # We do not increment i here.
                j = lps[j - 1]
            else:
                lps[i] = 0
                i += 1
    return lps


def kmp_search(text: str, pattern: str) -> list[int]:
    """
    Executes the KMP search algorithm to find all occurrences of a pattern.

    Args:
        text (str): The DNA sequence to search within.
        pattern (str): The motif/pattern to search for.

    Returns:
        list[int]: A list of starting indices where the pattern occurs.
                   Returns empty list if pattern is not found.
    """
    # 1. Edge Case: Empty pattern or pattern longer than text
    if not pattern:
        return []
    if len(pattern) > len(text):
        return []

    n = len(text)
    m = len(pattern)

    # 2. Preprocess
    lps = compute_lps_table(pattern)

    # 3. Search
    i = 0  # Index for text
    j = 0  # Index for pattern
    occurrences = []

    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1

        if j == m:
            # Match found!
            occurrences.append(i - j)
            # Prepare for next match (overlap handling)
            j = lps[j - 1]
        elif i < n and pattern[j] != text[i]:
            # Mismatch
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
                
    return occurrences

def count_comparisons_kmp(text: str, pattern: str) -> tuple[list[int], int]:
    """
    Executes KMP and counts the number of character comparisons performed.
    
    This is used for empirical analysis to demonstrate the efficiency
    difference between KMP (approx 2N comparisons) and Boyer-Moore (sub-linear).
    
    Args:
        text (str): The DNA sequence.
        pattern (str): The motif to search for.
        
    Returns:
        tuple: (list of matches, total_comparisons_made)
    """
    if not pattern or len(pattern) > len(text):
        return [], 0

    n = len(text)
    m = len(pattern)
    lps = compute_lps_table(pattern)

    i = 0
    j = 0
    occurrences = []
    comparisons = 0

    while i < n:
        # We are about to compare text[i] and pattern[j]
        comparisons += 1
        
        if pattern[j] == text[i]:
            i += 1
            j += 1
            
            if j == m:
                occurrences.append(i - j)
                j = lps[j - 1]
        else:
            # Mismatch happened (we already counted the comparison above)
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
                
    return occurrences, comparisons