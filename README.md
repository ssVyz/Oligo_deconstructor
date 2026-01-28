# Oligo Sequence Deconstruction Tool for qPCR Design

A professional Python/Tkinter application for analyzing and consolidating oligonucleotide sequences for qPCR primer design.

## Features

- **Parse FASTA files** with hundreds or thousands of sequences
- **Quality filtering**: Automatically removes sequences with ambiguities, gaps, or incorrect lengths
- **Multiple analysis modes**:
  - Find all unique sequence variants
  - Find minimum variants using up to N ambiguities (IUPAC codes)
  - Find top N most frequent variants
- **G-A wobble base pairing support**: Option to treat G as A (G-T mispairing not considered a mismatch)
- **Professional output**: Tab-separated results ready for Excel import
- **Copy-friendly**: One-click copy to clipboard

## Installation

### Prerequisites

- Python 3.8 or higher
- Tkinter (usually included with Python)

### Linux (Ubuntu/Debian)

```bash
# Install Tkinter if not already installed
sudo apt-get install python3-tk

# Clone or download the files
# No additional pip packages required
```

### Windows / macOS

Tkinter is typically included with Python installations from python.org.

### Running the Application

```bash
python oligo_analyzer.py
```

## Usage

### 1. Load Sequences

**Input Tab**: Load your sequences in FASTA format:

```
>Sequence1
AAAAAAAAA
>Sequence2
ATAAAAAAA
>Sequence3
ATAAAAAAA
```

You can:
- Click "Load FASTA File" to open a file
- Click "Paste from Clipboard" to paste sequences
- Click "Load Example" to see sample data

### 2. Configure Analysis

**Analysis Options Tab**: Select your analysis type and configure options:

#### Analysis Types

| Type | Description |
|------|-------------|
| **All unique variants** | Lists every unique sequence with frequency counts |
| **Minimum variants with N ambiguities** | Consolidates sequences using IUPAC ambiguity codes |
| **Top N variants** | Returns only the N most frequent variants |

#### Options

- **Maximum ambiguities**: Max number of ambiguous positions per variant (for min variants mode)
- **Number of top variants**: How many variants to return (for top N mode)
- **Treat G as A**: Enable G-T wobble base pairing tolerance

### 3. View Results

**Results Tab**: View, copy, or export your results.

Results are formatted with tabs for easy pasting into Excel:

```
Variant 1    AAA AAA AAA    3    50.0%
Variant 2    ATA AAA AAA    2    33.3%
Variant 3    AAA AGA AAA    1    16.7%
```

## IUPAC Ambiguity Codes

The tool uses standard IUPAC nucleotide codes:

| Code | Bases | Meaning |
|------|-------|---------|
| A | A | Adenine |
| C | C | Cytosine |
| G | G | Guanine |
| T | T | Thymine |
| R | A, G | Purine |
| Y | C, T | Pyrimidine |
| S | G, C | Strong |
| W | A, T | Weak |
| K | G, T | Keto |
| M | A, C | Amino |
| B | C, G, T | Not A |
| D | A, G, T | Not C |
| H | A, C, T | Not G |
| V | A, C, G | Not T |
| N | A, C, G, T | Any |

## Examples

### Example 1: All Variants Without Ambiguities

**Input:**
```
>Seq1: AAA AAA AAA
>Seq2: ATA AAA AAA
>Seq3: ATA AAA AAA
>Seq4: AAA AAA AAA
>Seq5: AAA AAA AAA
>Seq6: AAA AGA AAA
```

**Output:**
```
Variant 1    AAA AAA AAA    3    50.0%
Variant 2    ATA AAA AAA    2    33.3%
Variant 3    AAA AGA AAA    1    16.7%
```

### Example 2: Minimum Variants with 1 Ambiguity

**Result:**
```
Variant 1    AWA AAA AAA    5    83.3%
Variant 2    AAA AGA AAA    1    16.7%
```

(W = A or T, covering both AAA and ATA at position 2)

### Example 3: With G-A Wobble Enabled

When "Treat G as A" is enabled, G-T wobble base pairing is tolerated:

**Result:**
```
Variant 1    AWA AGA AAA    6    100%
```

## Quality Filtering

The tool automatically filters input sequences:

| Issue | Action |
|-------|--------|
| Ambiguous bases (R, Y, N, etc.) | Removed |
| Gap characters (-, .) | Removed |
| Wrong length (minorities) | Removed |

A detailed quality report shows how many sequences were removed and why.

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| Ctrl+O | Open FASTA file |
| Ctrl+R | Run analysis |

## Output Format

Results are tab-separated for easy import into spreadsheet applications:

```
Variant    Sequence    Count    Percentage
Variant 1    AAA AAA AAA    3    50.0%
```

Copy the results and paste directly into Excel - columns will be properly aligned.

## Extending the Tool

The code is designed for easy extension. To add new analysis options:

1. Add new UI controls in `create_analysis_tab()`
2. Implement the analysis method in `OligoAnalyzer` class
3. Add handling in `run_analysis()` method


## Troubleshooting

### "No module named 'tkinter'"

Install Tkinter for your system:
- Ubuntu/Debian: `sudo apt-get install python3-tk`
- Fedora: `sudo dnf install python3-tkinter`
- macOS with Homebrew: `brew install python-tk`

### Slow performance with many sequences

The minimum variants algorithm uses a greedy approach optimized for large datasets. For very large datasets (>10,000 sequences), consider:
- Using "Top N variants" mode first for quick overview
- Reducing the maximum ambiguities parameter

### Results don't paste correctly into Excel

Ensure you're using "Paste Special" > "Text" in Excel, or paste into a plain text editor first.
