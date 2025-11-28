"""
Boyer-Moore String Matching Algorithm.

Implementation of the exact pattern matching algorithm using 
Bad Character and Good Suffix heuristics.
"""

def bad_character_table(pattern: str) -> dict[str, int]:
    """
    Preprocesses the pattern to create the Bad Character heuristic table.
    
    Args:
        pattern (str): The search pattern.
    
    Returns:
        dict[str, int]: Map of character to its rightmost position (0-indexed).
    """
    bad_char = {}
    for i in range(len(pattern)):
        bad_char[pattern[i]] = i
    return bad_char


def good_suffix_table(pattern: str) -> list[int]:
    """
    Preprocesses the pattern to create the Good Suffix heuristic table.
    
    Args:
        pattern (str): The search pattern.
    
    Returns:
        list[int]: Array where shift[i] is the shift distance for mismatch at i.
    """
    m = len(pattern)
    shift = [m] * m
    border = [0] * (m + 1)
    
    # Case 2: Matched suffix matches a prefix of pattern
    i = m
    j = m + 1
    border[i] = j
    
    while i > 0:
        while j <= m and pattern[i - 1] != pattern[j - 1]:
            if shift[j - 1] == m:
                shift[j - 1] = j - i
            j = border[j]
        i -= 1
        j -= 1
        border[i] = j
    
    # Case 1: Matched suffix occurs elsewhere in pattern
    j = border[0]
    for i in range(m):
        if shift[i] == m:
            shift[i] = j
        if i == j:
            j = border[j]
            
    return shift


def boyer_moore_search(text: str, pattern: str) -> list[int]:
    """
    Executes the Boyer-Moore search algorithm.
    
    Args:
        text (str): The DNA sequence to search within.
        pattern (str): The motif to find.
    
    Returns:
        list[int]: List of starting indices where pattern occurs.
    """
    if not pattern:
        return []
    
    n = len(text)
    m = len(pattern)
    
    if m > n:
        return []
    
    # Preprocessing
    bad_char = bad_character_table(pattern)
    good_suffix = good_suffix_table(pattern)
    
    matches = []
    shift_amount = 0
    
    while shift_amount <= n - m:
        j = m - 1
        
        # Compare from right to left
        while j >= 0 and pattern[j] == text[shift_amount + j]:
            j -= 1
        
        if j < 0:
            # Match found
            matches.append(shift_amount)
            # Shift pattern
            shift_amount += good_suffix[0] if shift_amount + m < n else 1
        else:
            # Mismatch
            bad_char_shift = j - bad_char.get(text[shift_amount + j], -1)
            good_suffix_shift = good_suffix[j]
            shift_amount += max(bad_char_shift, good_suffix_shift)
            
    return matches


def count_comparisons_boyer_moore(text: str, pattern: str) -> tuple[list[int], int]:
    """
    Runs Boyer-Moore and counts character comparisons.
    
    Used for empirical analysis in the project report to demonstrate
    sub-linear behavior compared to KMP or Naive approaches.
    
    Returns:
        tuple: (list of matches, total_comparisons_made)
    """
    if not pattern:
        return [], 0
        
    n = len(text)
    m = len(pattern)
    
    if m > n:
        return [], 0
    
    bad_char = bad_character_table(pattern)
    good_suffix = good_suffix_table(pattern)
    
    matches = []
    comparisons = 0
    shift_amount = 0
    
    while shift_amount <= n - m:
        j = m - 1
        
        while j >= 0:
            comparisons += 1
            if pattern[j] != text[shift_amount + j]:
                break
            j -= 1
        
        if j < 0:
            matches.append(shift_amount)
            shift_amount += good_suffix[0] if shift_amount + m < n else 1
        else:
            bad_char_shift = j - bad_char.get(text[shift_amount + j], -1)
            good_suffix_shift = good_suffix[j]
            shift_amount += max(bad_char_shift, good_suffix_shift)
            
    return matches, comparisons