"""
Unit Tests for Boyer-Moore Algorithm
Tests correctness, edge cases, and DNA-specific scenarios

Author: Team Subproblem Solvers
"""

import sys
import os
import unittest

# Add 'src' directory to path so we can import modules
sys.path.append(os.path.join(os.getcwd(), 'src'))

from algorithms.boyer_moore import (
    boyer_moore_search,
    bad_character_table,
    good_suffix_table,
    count_comparisons_boyer_moore
)

class TestBadCharacterTable(unittest.TestCase):
    """Test the bad character heuristic preprocessing."""
    
    def test_simple_pattern(self):
        """Test bad character table for simple pattern."""
        result = bad_character_table("GATC")
        self.assertEqual(result, {'G': 0, 'A': 1, 'T': 2, 'C': 3})
    
    def test_repeated_characters(self):
        """Test pattern with repeated characters (rightmost should be stored)."""
        result = bad_character_table("GCAGAGAG")
        self.assertEqual(result['G'], 7)
        self.assertEqual(result['A'], 6)
        self.assertEqual(result['C'], 1)
    
    def test_single_character(self):
        """Test pattern with single character."""
        result = bad_character_table("A")
        self.assertEqual(result, {'A': 0})
    
    def test_dna_pattern(self):
        """Test with typical DNA motif."""
        result = bad_character_table("ATCG")
        self.assertEqual(len(result), 4)


class TestGoodSuffixTable(unittest.TestCase):
    """Test the good suffix heuristic preprocessing."""
    
    def test_table_length(self):
        pattern = "GATC"
        result = good_suffix_table(pattern)
        self.assertEqual(len(result), len(pattern))
    
    def test_short_pattern(self):
        pattern = "AT"
        result = good_suffix_table(pattern)
        self.assertEqual(len(result), 2)


class TestBoyerMooreSearch(unittest.TestCase):
    """Test the main Boyer-Moore search algorithm."""
    
    def test_simple_match(self):
        """Test finding a single occurrence."""
        text = "ATCGATCG"
        pattern = "GATC"
        # Correct index is 3 (A=0, T=1, C=2, G=3)
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [3])
    
    def test_multiple_matches(self):
        """Test finding multiple occurrences."""
        text = "ATCGATCGATCG"
        pattern = "GATC"
        # Matches at index 3 and index 7
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [3, 7])
    
    def test_no_match(self):
        """Test when pattern doesn't occur."""
        text = "ATCGATCG"
        pattern = "GGGG"
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [])
    
    def test_pattern_longer_than_text(self):
        """Test when pattern is longer than text."""
        text = "ATCG"
        pattern = "ATCGATCG"
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [])
    
    def test_empty_pattern(self):
        """Test with empty pattern."""
        text = "ATCG"
        pattern = ""
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [])
    
    def test_pattern_at_start(self):
        """Test pattern at the beginning of text."""
        text = "GATCATCG"
        pattern = "GATC"
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [0])
    
    def test_pattern_at_end(self):
        """Test pattern at the end of text."""
        text = "ATCGGATC"
        pattern = "GATC"
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [4])
    
    def test_overlapping_patterns(self):
        """Test with overlapping pattern occurrences."""
        text = "AAAAAAA"
        pattern = "AAA"
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [0, 1, 2, 3, 4])
    
    def test_entire_text_is_pattern(self):
        """Test when the entire text is the pattern."""
        text = "ATCG"
        pattern = "ATCG"
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [0])
    
    def test_repeated_pattern(self):
        """Test text with repeated pattern."""
        text = "GATCGATCGATC"
        pattern = "GATC"
        # GATC starts at 0, 4, 8
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [0, 4, 8])


class TestDNASpecificCases(unittest.TestCase):
    """Test DNA-specific scenarios."""
    
    def test_dna_motif_tataa_box(self):
        """Test finding TATA box."""
        text = "CGCGTATAACGCGTATAAGGC"
        pattern = "TATAA"
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [4, 13])
    
    def test_start_codon(self):
        """Test finding start codon ATG."""
        text = "CGCGATGCCGATGAAATG"
        pattern = "ATG"
        # CGCG(4)ATG(7)CCG(10)ATG(13)AA(15)ATG
        result = boyer_moore_search(text, pattern)
        self.assertEqual(result, [4, 10, 15])
    
    def test_gc_rich_region(self):
        """Test pattern in GC-rich region."""
        text = "GCGCGCGCGCGC"
        pattern = "GCGC"
        result = boyer_moore_search(text, pattern)
        # Should find multiple overlapping matches
        self.assertTrue(len(result) > 0)
        self.assertIn(0, result)


class TestComparisonCounting(unittest.TestCase):
    """Test the comparison counting variant."""
    
    def test_returns_matches_and_count(self):
        """Test that function returns both matches and comparison count."""
        text = "ATCGATCG"
        pattern = "GATC"
        matches, comparisons = count_comparisons_boyer_moore(text, pattern)
        
        self.assertIsInstance(matches, list)
        self.assertIsInstance(comparisons, int)
        self.assertEqual(matches, [3])
        self.assertGreater(comparisons, 0)
    
    def test_comparison_count_reasonable(self):
        """Test that comparison count is reasonable."""
        text = "A" * 1000
        pattern = "AAA"
        matches, comparisons = count_comparisons_boyer_moore(text, pattern)
        self.assertGreater(comparisons, 0)

def run_tests():
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])
    runner.run(suite)

if __name__ == "__main__":
    run_tests()