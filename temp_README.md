# Genomic Pattern Matching Project

## Team Subproblem Solvers
A comparative analysis of string matching algorithms on Genomic Data.

## Project Structure
*   `src/algorithms/`: Contains the algorithm logic.
*   `src/utils/`: Contains data generators.
*   `Analysis.ipynb`: Main benchmarking and visualization notebook.
*   `test_*.py`: Unit tests for ensuring correctness.


## How to run the analysis
1.  Open `Analysis.ipynb` in VS Code or Jupyter.
2.  Click "Run All".
3.  Graphs for Scalability and Efficiency will be generated automatically.

## For Team Members
**To add your algorithm (e.g., Suffix Tree):**
1.  Create your file in `src/algorithms/suffix_tree.py`.
2.  Implement your search function.
3.  Open `Analysis.ipynb`.
4.  In **Cell 2**, import your function and add it to the `algorithms` dictionary.
5.  Run the notebook to see your algorithm compared against KMP and Boyer-Moore!