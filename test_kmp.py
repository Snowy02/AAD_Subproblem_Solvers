"""
Unit Tests for KMP Algorithm
Tests correctness, edge cases, and DNA-specific scenarios.

Author: Team Subproblem Solvers
"""

import sys
import os
import unittest

# Add 'src' directory to path so we can import modules
sys.path.append(os.path.join(os.getcwd(), 'src'))

from algorithms.kmp import kmp_search, compute_lps_table, count_comparisons_kmp

class TestLPSTable(unittest.TestCase):
    """Test the Longest Prefix-Suffix (LPS) preprocessing."""

    def test_basic_pattern(self):
        """Test LPS for simple repetitive pattern."""
        self.assertEqual(compute_lps_table("AAAA"), [0, 1, 2, 3])

    def test_no_repetition(self):
        """Test pattern with no repeating characters."""
        self.assertEqual(compute_lps_table("ABCDE"), [0, 0, 0, 0, 0])

    def test_complex_dna(self):
        """Test standard DNA mix."""
        # Pattern: AABAACAABAA
        expected = [0, 1, 0, 1, 2, 0, 1, 2, 3, 4, 5]
        self.assertEqual(compute_lps_table("AABAACAABAA"), expected)

    def test_single_char(self):
        self.assertEqual(compute_lps_table("A"), [0])


class TestKMPSearch(unittest.TestCase):
    """Test the main KMP search algorithm."""
    
    def test_simple_match(self):
        """Test finding a single occurrence."""
        text = "ATCGATCG"
        pattern = "GATC"
        # Correct index is 3 (A=0, T=1, C=2, G=3)
        result = kmp_search(text, pattern)
        self.assertEqual(result, [3])
    
    def test_multiple_matches(self):
        """Test finding multiple occurrences."""
        text = "ATCGATCGATCG"
        pattern = "GATC"
        # Matches at index 3 and index 7
        result = kmp_search(text, pattern)
        self.assertEqual(result, [3, 7])
    
    def test_no_match(self):
        """Test when pattern doesn't occur."""
        text = "ATCGATCG"
        pattern = "GGGG"
        result = kmp_search(text, pattern)
        self.assertEqual(result, [])
    
    def test_pattern_longer_than_text(self):
        """Test when pattern is longer than text."""
        text = "ATCG"
        pattern = "ATCGATCG"
        result = kmp_search(text, pattern)
        self.assertEqual(result, [])
    
    def test_empty_pattern(self):
        """Test with empty pattern (Should not crash)."""
        text = "ATCG"
        pattern = ""
        result = kmp_search(text, pattern)
        self.assertEqual(result, [])
    
    def test_overlapping_patterns(self):
        """Test with overlapping pattern occurrences (Crucial for KMP)."""
        text = "AAAAAAA"
        pattern = "AAA"
        # Indices: 012, 123, 234, 345, 456
        result = kmp_search(text, pattern)
        self.assertEqual(result, [0, 1, 2, 3, 4])


class TestDNASpecificCases(unittest.TestCase):
    """Test DNA-specific scenarios relevant to the project."""
    
    def test_dna_motif_tataa_box(self):
        """Test finding TATA box (common promoter sequence)."""
        text = "CGCGTATAACGCGTATAAGGC"
        pattern = "TATAA"
        result = kmp_search(text, pattern)
        self.assertEqual(result, [4, 13])
    
    def test_start_codon(self):
        """Test finding start codon ATG."""
        text = "CGCGATGCCGATGAAATG"
        pattern = "ATG"
        # 0123 456 789 012 34 567
        # CGCG ATG CCG ATG AA ATG
        result = kmp_search(text, pattern)
        self.assertEqual(result, [4, 10, 15])
    
    def test_gc_rich_region(self):
        """Test pattern in GC-rich region."""
        text = "GCGCGCGCGCGC"
        pattern = "GCGC"
        result = kmp_search(text, pattern)
        # Should find overlapping matches
        self.assertTrue(len(result) > 0)
        self.assertIn(0, result)
        self.assertIn(2, result)
class TestKMPComparisonCounting(unittest.TestCase):
    """Test the comparison counting variant."""
    
    def test_counts_returned(self):
        text = "ATCGATCG"
        pattern = "GATC"
        matches, comps = count_comparisons_kmp(text, pattern)
        self.assertEqual(matches, [3])
        self.assertGreater(comps, 0)
        
    def test_efficiency(self):
        # KMP should be roughly linear. 
        # For 1000 chars, comparisons should be close to 1000-2000, not 1,000,000.
        text = "A" * 1000
        pattern = "AAA"
        matches, comps = count_comparisons_kmp(text, pattern)
        self.assertLess(comps, 3000)
def run_tests():
    """Run KMP test suite with clean formatted output."""
    print("\n=== Running KMP Algorithm Tests ===")
    runner = unittest.TextTestRunner(verbosity=2)
    suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])
    result = runner.run(suite)
    print("=== KMP Testing Complete ===\n")

if __name__ == "__main__":
    run_tests()
