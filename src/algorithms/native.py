"""
Native Python String Search Benchmarking Wrapper.

This module wraps Python's built-in 'find' method and 're' module
to serve as a baseline for performance comparisons.
"""
import re

def native_search(text: str, pattern: str) -> list[int]:
    """
    Finds all occurrences using Python's highly optimized C-based .find() method.
    Acts as a 'Baseline' for Exact Matching speed.
    """
    matches = []
    start = 0
    while True:
        idx = text.find(pattern, start)
        if idx == -1:
            break
        matches.append(idx)
        start = idx + 1
    return matches

def regex_search(text: str, pattern: str) -> list[int]:
    """
    Finds all occurrences using Python's 're' module.
    Acts as a 'Baseline' for general pattern matching.
    """
    # Escape pattern to treat it as a literal string for fair comparison
    # (unless specifically testing regex features)
    escaped_pattern = re.escape(pattern)
    return [m.start() for m in re.finditer(escaped_pattern, text)]