"""
Levenshtein Distance Search (Edit Distance)
Approximate string matching for genomic analysis.

1. Standard Levenshtein distance (Wagner–Fischer DP)
2. Naive sliding-window approximate search
3. Optimized rolling-DP search
4. Diagonal/banded DP search (Ukkonen-style)
5. Comparison-count version for empirical evaluation

"""


def levenshtein_distance(a: str, b: str) -> int:

    n, m = len(a), len(b)

    # Edge cases
    if n == 0:
        return m
    if m == 0:
        return n

    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        dp[i][0] = i
    for j in range(m + 1):
        dp[0][j] = j

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
 
    n, m = len(text), len(pattern)

    if m == 0:
        return list(range(n + 1))  
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

    n, m = len(text), len(pattern)

    if m == 0:
        return list(range(n + 1))
    if n < m:
        return []

    matches = []

    for start in range(n - m + 1):
        substring = text[start:start + m]

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
  
    n, m = len(text), len(pattern)

    if m == 0:
        return list(range(n + 1))
    if n < m:
        return []

    matches = []

    prev = list(range(m + 1))  

    for i in range(1, n + 1):
        curr = [0] * (m + 1)

       
        start_j = max(1, i - max_distance)
        end_j = min(m, i + max_distance)

        curr[0] = i 

        for j in range(start_j, end_j + 1):
            cost = 0 if text[i - 1] == pattern[j - 1] else 1

            ins = curr[j - 1] + 1
            delete = prev[j] + 1
            subst = prev[j - 1] + cost

            curr[j] = min(ins, delete, subst)

        for j in range(1, start_j):
            curr[j] = max_distance + 5
        for j in range(end_j + 1, m + 1):
            curr[j] = max_distance + 5

        if i >= m and curr[m] <= max_distance:
            matches.append(i - m)

        prev = curr

    return matches



def count_comparisons_levenshtein(text: str, pattern: str, max_distance: int):
   
    n, m = len(text), len(pattern)
    comparisons = 0

    if m == 0:
        return list(range(n + 1)), 0
    if n < m:
        return [], 0

    matches = []

    for start in range(n - m + 1):
        prev = list(range(m + 1))
        curr = [0] * (m + 1)

        for i in range(1, m + 1):
            curr[0] = i
            for j in range(1, m + 1):
                comparisons += 1
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
