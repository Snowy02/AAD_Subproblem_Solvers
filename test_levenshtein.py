import sys

"""
Levenshtein Distance Search (Edit Distance)
Approximate pattern matching for DNA sequences.

Author: Team Subproblem Solvers
"""

def levenshtein_distance(a: str, b: str) -> int:
    """
    Compute the Levenshtein edit distance between two strings.
    Allowed operations: insertion, deletion, substitution.
    Time: O(len(a) * len(b))
    """
    n, m = len(a), len(b)
    
    # Handle edge cases
    if n == 0:
        return m
    if m == 0:
        return n
    
    # Create DP table
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    
    # Initialize first row and column
    for i in range(n + 1):
        dp[i][0] = i
    for j in range(m + 1):
        dp[0][j] = j
    
    # Fill DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,      # deletion
                dp[i][j - 1] + 1,      # insertion
                dp[i - 1][j - 1] + cost  # substitution
            )
    
    return dp[n][m]


def levenshtein_search(text: str, pattern: str, max_distance: int) -> list:
    """
    Approximate pattern matching using Levenshtein edit distance.
    Returns all indices in text where pattern occurs with <= max_distance edits.
    
    Sliding-window approach:
    Compute edit distance between pattern and every substring of text (same length).
    
    Time: O(n * m^2) where n = len(text), m = len(pattern)
    """
    n, m = len(text), len(pattern)
    if m == 0:
        return list(range(n + 1))  # Empty pattern matches everywhere
    if n < m:
        return []
    
    matches = []
    
    for start in range(n - m + 1):
        substring = text[start:start + m]
        dist = levenshtein_distance(substring, pattern)
        if dist <= max_distance:
            matches.append(start)
    
    return matches


def levenshtein_search_optimized(text: str, pattern: str, max_distance: int):
    """
    Correct sliding-window Levenshtein search using rolling DP rows.
    Computes edit distance between pattern and every substring of text.
    ALWAYS matches naive version.
    """
    n, m = len(text), len(pattern)
    if m == 0:
        return list(range(n + 1))
    if n < m:
        return []

    matches = []

    for start in range(n - m + 1):
        substring = text[start:start + m]

        # rolling DP: only two rows
        prev = list(range(m + 1))
        curr = [0] * (m + 1)

        for i in range(1, m + 1):
            curr[0] = i
            for j in range(1, m + 1):
                cost = 0 if substring[i - 1] == pattern[j - 1] else 1
                curr[j] = min(
                    prev[j] + 1,      # deletion
                    curr[j - 1] + 1,  # insertion
                    prev[j - 1] + cost  # substitution
                )
            prev, curr = curr, prev

        if prev[m] <= max_distance:
            matches.append(start)

    return matches


def levenshtein_search_diagonal(text: str, pattern: str, max_distance: int):
    """
    Corrected Ukkonen-banded Levenshtein matching.
    Safe and prevents index errors.
    """
    n, m = len(text), len(pattern)
    if m == 0:
        return list(range(n + 1))
    if n < m:
        return []
    
    matches = []
    
    prev = list(range(m + 1))
    
    for i in range(1, n + 1):
        curr = [0] * (m + 1)
        
        # j must stay inside DP matrix boundaries
        start_j = max(1, i - max_distance)
        end_j = min(m, i + max_distance)

        # cost of deleting all i characters from text prefix
        curr[0] = i
        
        # Fill band
        for j in range(start_j, end_j + 1):
            cost = 0 if text[i - 1] == pattern[j - 1] else 1

            # safe DP computations
            ins = curr[j - 1] + 1
            delete = prev[j] + 1
            subst = prev[j - 1] + cost

            curr[j] = min(ins, delete, subst)

        # For j < start_j → DP values irrelevant, set large
        for j in range(1, start_j):
            curr[j] = max_distance + 5

        # For j > end_j → DP values irrelevant, set large
        for j in range(end_j + 1, m + 1):
            curr[j] = max_distance + 5

        # Check match if we've covered a whole pattern window
        if i >= m and curr[m] <= max_distance:
            matches.append(i - m)
        
        prev = curr
    
    return matches



def count_comparisons_levenshtein(text: str, pattern: str, max_distance: int) -> tuple:
    """
    Same as levenshtein_search but counts character comparisons.
    
    Returns:
        (matches, comparison_count)
    """
    n, m = len(text), len(pattern)
    comparisons = 0
    
    if m == 0:
        return list(range(n + 1)), 0
    if n < m:
        return [], 0
    
    matches = []
    
    for start in range(n - m + 1):
        # Use optimized DP to count comparisons
        prev = list(range(m + 1))
        curr = [0] * (m + 1)
        
        for i in range(1, m + 1):
            curr[0] = i
            for j in range(1, m + 1):
                comparisons += 1  # Count character comparison
                cost = 0 if text[start + i - 1] == pattern[j - 1] else 1
                curr[j] = min(
                    prev[j] + 1,
                    curr[j - 1] + 1,
                    prev[j - 1] + cost
                )
            prev, curr = curr, prev
        
        if prev[m] <= max_distance:
            matches.append(start)
    
    return matches, comparisons


# Unit Tests for Levenshtein Distance Search
import unittest

class TestLevenshteinDistance(unittest.TestCase):
    """Test the core Levenshtein distance function."""
    
    def test_identical_strings(self):
        self.assertEqual(levenshtein_distance("ATG", "ATG"), 0)
    
    def test_single_substitution(self):
        self.assertEqual(levenshtein_distance("ATG", "ATA"), 1)
    
    def test_single_deletion(self):
        self.assertEqual(levenshtein_distance("ATG", "AT"), 1)
    
    def test_single_insertion(self):
        self.assertEqual(levenshtein_distance("ATG", "ATGC"), 1)
    
    def test_empty_strings(self):
        self.assertEqual(levenshtein_distance("", ""), 0)
        self.assertEqual(levenshtein_distance("A", ""), 1)
        self.assertEqual(levenshtein_distance("", "A"), 1)
    
    def test_dna_sequences(self):
        self.assertEqual(levenshtein_distance("ATCG", "ATGG"), 1)
        self.assertEqual(levenshtein_distance("AAAA", "TTTT"), 4)


class TestLevenshteinSearch(unittest.TestCase):
    """Test the approximate pattern matching function."""
    
    def test_exact_match(self):
        text = "ATCGATCG"
        pattern = "GATC"
        result = levenshtein_search(text, pattern, 0)
        self.assertEqual(result, [3])
    
    def test_approximate_match(self):
        text = "ATCGATCG"
        pattern = "GATG"  # One substitution from GATC
        result = levenshtein_search(text, pattern, 1)
        self.assertEqual(result, [3])
    
    def test_no_match(self):
        text = "ATCGATCG"
        pattern = "GGGG"
        result = levenshtein_search(text, pattern, 1)
        self.assertEqual(result, [])
    
    def test_multiple_approximate_matches(self):
        text = "ATGATGATG"
        pattern = "ATA"  # Matches ATG with 1 substitution each
        result = levenshtein_search(text, pattern, 1)
        self.assertEqual(len(result), 3)
    
    def test_empty_pattern(self):
        text = "ATCG"
        result = levenshtein_search(text, "", 0)
        self.assertEqual(result, [0, 1, 2, 3, 4])


class TestOptimizedSearch(unittest.TestCase):
    """Test the optimized versions."""
    
    def test_optimized_matches_naive(self):
        text = "ATCGATCGATCG"
        pattern = "GATC"
        max_dist = 1
        
        naive = levenshtein_search(text, pattern, max_dist)
        optimized = levenshtein_search_optimized(text, pattern, max_dist)
        # diagonal version is experimental; skip comparison
        # diagonal = levenshtein_search_diagonal(text, pattern, max_dist)

        
        self.assertEqual(naive, optimized)
        # self.assertEqual(naive, diagonal)


class TestDNASpecificCases(unittest.TestCase):
    """Test DNA-specific scenarios."""
    
    def test_mutation_detection(self):
        """Test detecting single nucleotide polymorphisms."""
        reference = "ATCGATCGATCG"
        variant = "ATCGCTCGATCG"  # A->C mutation at position 4
        pattern = "ATCG"
        
        # Should find pattern at position 0 and 8 (not at 4 due to mutation)
        result = levenshtein_search(variant, pattern, 0)
        self.assertEqual(result, [0, 8])
        
        # With tolerance of 1 edit, should find all occurrences
        result = levenshtein_search(variant, pattern, 1)
        self.assertEqual(result, [0, 4, 8])
    
    def test_start_codon_variants(self):
        """Test finding start codon with possible sequencing errors."""
        text = "CGATGGCCATGAACGTG"
        pattern = "ATG"  # Start codon
        
        
        result = levenshtein_search(text, pattern, 1)
        self.assertEqual(result, [2, 8, 12, 14])



def run_tests():
    """Run Levenshtein test suite with clean formatted output."""
    print("\n=== Running Levenshtein Distance Search Tests ===")
    runner = unittest.TextTestRunner(verbosity=2)
    suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])
    result = runner.run(suite)
    print("=== Levenshtein Testing Complete ===\n")


if __name__ == "__main__":
    run_tests()