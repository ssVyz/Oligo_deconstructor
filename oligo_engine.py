#!/usr/bin/env python3
"""
Oligo Sequence Analysis Engine
==============================
Core analysis functionality without GUI dependencies.

Can be imported independently for scripting or integrated with GUI.
"""

from collections import Counter
from typing import Dict, List, Tuple, Set, Optional, Callable
from dataclasses import dataclass, field


# ============================================================================
# IUPAC AMBIGUITY CODES
# ============================================================================
IUPAC_CODES = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

# Reverse lookup: bases -> IUPAC code
BASES_TO_IUPAC = {frozenset(v): k for k, v in IUPAC_CODES.items()}

# Complement mapping for reverse complement
COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'R': 'Y', 'Y': 'R',  # R(A/G) <-> Y(C/T)
    'S': 'S', 'W': 'W',  # S(G/C) and W(A/T) are self-complementary
    'K': 'M', 'M': 'K',  # K(G/T) <-> M(A/C)
    'B': 'V', 'V': 'B',  # B(C/G/T) <-> V(A/C/G)
    'D': 'H', 'H': 'D',  # D(A/G/T) <-> H(A/C/T)
    'N': 'N',            # N is self-complementary
}

AMBIGUOUS_BASES = set('RYSWKMBDHVN')
VALID_BASES = set('ACGTRYSWKMBDHVN')
GAP_CHARS = set('-.')


# ============================================================================
# DATA CLASSES
# ============================================================================
@dataclass
class AnalysisResult:
    """Container for analysis results"""
    variants: List[Tuple[str, int, float]] = field(default_factory=list)
    total_sequences: int = 0
    coverage: float = 0.0
    uncovered_count: int = 0
    uncovered_percentage: float = 0.0
    ambiguity_count: int = 0
    message: str = ""


@dataclass
class QualityReport:
    """Container for quality filtering report"""
    original_count: int = 0
    valid_count: int = 0
    removed_ambiguous: int = 0
    removed_gaps: int = 0
    removed_wrong_length: int = 0
    majority_length: int = 0
    valid_sequences: List[str] = field(default_factory=list)


# ============================================================================
# SEQUENCE ANALYSIS ENGINE
# ============================================================================
class OligoAnalyzer:
    """Core analysis engine for oligo sequence processing"""
    
    def __init__(self):
        self.sequences: List[str] = []
        self.quality_report: Optional[QualityReport] = None
    
    def parse_fasta(self, text: str) -> List[str]:
        """Parse FASTA format text and extract sequences"""
        sequences = []
        current_seq = []
        
        lines = text.strip().split('\n')
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq).upper().replace(' ', ''))
                    current_seq = []
            else:
                # Remove spaces and convert to uppercase
                current_seq.append(line.upper().replace(' ', ''))
        
        # Don't forget the last sequence
        if current_seq:
            sequences.append(''.join(current_seq).upper().replace(' ', ''))
        
        # If no FASTA headers found, treat each line as a sequence
        if not sequences and lines:
            for line in lines:
                line = line.strip().upper().replace(' ', '')
                if line and not line.startswith('>'):
                    sequences.append(line)
        
        return sequences
    
    def quality_filter(self, sequences: List[str]) -> QualityReport:
        """Filter sequences for quality issues"""
        report = QualityReport()
        report.original_count = len(sequences)
        
        if not sequences:
            return report
        
        # Determine majority length
        length_counts = Counter(len(seq) for seq in sequences)
        report.majority_length = length_counts.most_common(1)[0][0]
        
        valid = []
        for seq in sequences:
            # Check for gaps
            if any(c in GAP_CHARS for c in seq):
                report.removed_gaps += 1
                continue
            
            # Check for ambiguous bases
            if any(c in AMBIGUOUS_BASES for c in seq):
                report.removed_ambiguous += 1
                continue
            
            # Check for invalid characters
            if not all(c in VALID_BASES for c in seq):
                report.removed_ambiguous += 1
                continue
            
            # Check length
            if len(seq) != report.majority_length:
                report.removed_wrong_length += 1
                continue
            
            valid.append(seq)
        
        report.valid_sequences = valid
        report.valid_count = len(valid)
        return report
    
    def load_sequences(self, text: str) -> QualityReport:
        """Load and filter sequences from text"""
        raw_sequences = self.parse_fasta(text)
        self.quality_report = self.quality_filter(raw_sequences)
        self.sequences = self.quality_report.valid_sequences
        return self.quality_report
    
    def count_variants(self) -> AnalysisResult:
        """Count all unique sequence variants"""
        result = AnalysisResult()
        result.total_sequences = len(self.sequences)
        
        if not self.sequences:
            result.message = "No sequences to analyze"
            return result
        
        counts = Counter(self.sequences)
        total = len(self.sequences)
        
        # Sort by frequency (most common first)
        sorted_variants = counts.most_common()
        
        for seq, count in sorted_variants:
            percentage = (count / total) * 100
            result.variants.append((seq, count, percentage))
        
        result.coverage = 100.0
        return result
    
    def _bases_match(self, base1: str, base2: str, treat_g_as_a: bool = False) -> bool:
        """Check if two bases match (considering wobble pairing)"""
        if base1 == base2:
            return True
        if treat_g_as_a:
            # G-A wobble: treat G as matching A
            if (base1 == 'G' and base2 == 'A') or (base1 == 'A' and base2 == 'G'):
                return True
        return False
    
    def _get_ambiguity_code(self, bases: Set[str], treat_g_as_a: bool = False) -> str:
        """Get IUPAC ambiguity code for a set of bases"""
        if treat_g_as_a:
            # If treating G as A, we can use A to represent both
            if bases == {'A', 'G'} or bases == {'A'} or bases == {'G'}:
                return 'A'  # Use A since G-T wobble pairs with A
        
        frozen = frozenset(bases)
        return BASES_TO_IUPAC.get(frozen, 'N')
    
    def _count_ambiguities(self, seq: str) -> int:
        """Count number of ambiguous bases in a sequence"""
        return sum(1 for c in seq if c in AMBIGUOUS_BASES)
    
    def _sequence_matches_consensus(self, seq: str, consensus: str, treat_g_as_a: bool = False) -> bool:
        """Check if a sequence matches a consensus (with ambiguities)"""
        if len(seq) != len(consensus):
            return False
        
        for s, c in zip(seq, consensus):
            if c in IUPAC_CODES:
                allowed = IUPAC_CODES[c]
                if treat_g_as_a:
                    # Expand allowed bases for G-A wobble
                    if 'A' in allowed:
                        allowed = allowed | {'G'}
                    if 'G' in allowed:
                        allowed = allowed | {'A'}
                if s not in allowed:
                    return False
            else:
                if not self._bases_match(s, c, treat_g_as_a):
                    return False
        return True
    
    def _create_consensus(self, sequences: List[str], treat_g_as_a: bool = False,
                          exclude_n: bool = False) -> Tuple[str, int, bool]:
        """Create a consensus sequence with ambiguity codes.

        Returns:
            Tuple of (consensus_sequence, ambiguity_count, is_valid)
            is_valid is False if exclude_n=True and an N would be required
        """
        if not sequences:
            return "", 0, True

        seq_len = len(sequences[0])
        consensus = []
        ambiguity_count = 0

        for pos in range(seq_len):
            bases_at_pos = set(seq[pos] for seq in sequences)

            if treat_g_as_a:
                # If we have both A and G, treat them as equivalent
                if bases_at_pos == {'A', 'G'} or bases_at_pos == {'A'} or bases_at_pos == {'G'}:
                    consensus.append('A')
                    continue
                # If we have A or G mixed with others
                if 'A' in bases_at_pos and 'G' in bases_at_pos:
                    bases_at_pos = (bases_at_pos - {'G'}) | {'A'}

            if len(bases_at_pos) == 1:
                consensus.append(bases_at_pos.pop())
            else:
                code = self._get_ambiguity_code(bases_at_pos, treat_g_as_a)
                # Check if N is required but excluded
                if exclude_n and code == 'N':
                    return ''.join(consensus), ambiguity_count, False
                consensus.append(code)
                ambiguity_count += 1

        return ''.join(consensus), ambiguity_count, True
    
    def find_minimum_variants_greedy(self, max_ambiguities: int, treat_g_as_a: bool = False,
                                      exclude_n: bool = False,
                                      progress_callback: Optional[Callable[[str], None]] = None) -> AnalysisResult:
        """
        Find minimum set of variants using at most n ambiguities per variant.
        Uses a greedy set-cover approach for efficiency with large datasets.

        Args:
            max_ambiguities: Maximum number of ambiguous positions allowed per variant
            treat_g_as_a: If True, treat G and A as equivalent (G-A wobble)
            exclude_n: If True, do not allow N (any base) as an ambiguity code
            progress_callback: Optional callback for progress updates
        """
        result = AnalysisResult()
        result.total_sequences = len(self.sequences)

        if not self.sequences:
            result.message = "No sequences to analyze"
            return result

        # Get unique sequences with counts
        seq_counts = Counter(self.sequences)
        unique_seqs = list(seq_counts.keys())
        total = len(self.sequences)

        # Track which sequences are still uncovered
        uncovered = set(unique_seqs)
        variants = []

        iteration = 0
        while uncovered:
            iteration += 1
            if progress_callback:
                progress_callback(f"Finding variant {len(variants) + 1}... ({len(uncovered)} sequences remaining)")

            best_consensus = None
            best_coverage = set()
            best_ambiguities = float('inf')

            # Try to find the best consensus that covers the most uncovered sequences
            # Start with the most frequent uncovered sequence
            uncovered_counts = [(seq, seq_counts[seq]) for seq in uncovered]
            uncovered_counts.sort(key=lambda x: -x[1])

            # For each potential seed sequence
            for seed_seq, _ in uncovered_counts[:min(50, len(uncovered_counts))]:  # Limit seeds for performance
                # Find all sequences that could potentially be combined with this seed
                potential_group = [seed_seq]

                for other_seq in uncovered:
                    if other_seq == seed_seq:
                        continue

                    # Check if combining would exceed ambiguity limit or require N
                    combined = potential_group + [other_seq]
                    consensus, amb_count, is_valid = self._create_consensus(combined, treat_g_as_a, exclude_n)

                    if is_valid and amb_count <= max_ambiguities:
                        potential_group.append(other_seq)

                # Create consensus for this group
                consensus, amb_count, is_valid = self._create_consensus(potential_group, treat_g_as_a, exclude_n)

                # If consensus is invalid (requires N but excluded), skip this group
                if not is_valid:
                    continue

                # Check actual coverage (may include sequences not in potential_group)
                coverage = set()
                for seq in uncovered:
                    if self._sequence_matches_consensus(seq, consensus, treat_g_as_a):
                        coverage.add(seq)

                # Prefer larger coverage, then fewer ambiguities
                coverage_score = sum(seq_counts[s] for s in coverage)
                if coverage_score > sum(seq_counts[s] for s in best_coverage) or \
                   (coverage_score == sum(seq_counts[s] for s in best_coverage) and amb_count < best_ambiguities):
                    best_consensus = consensus
                    best_coverage = coverage
                    best_ambiguities = amb_count

            if best_consensus is None or not best_coverage:
                # Fallback: use the most frequent uncovered sequence as-is
                most_freq = max(uncovered, key=lambda x: seq_counts[x])
                best_consensus = most_freq
                best_coverage = {most_freq}
                best_ambiguities = 0

            # Record this variant
            count = sum(seq_counts[s] for s in best_coverage)
            percentage = (count / total) * 100
            variants.append((best_consensus, count, percentage))
            result.ambiguity_count = max(result.ambiguity_count, best_ambiguities)

            # Remove covered sequences
            uncovered -= best_coverage

        result.variants = variants
        result.coverage = 100.0
        return result
    
    def find_top_n_variants(self, n: int, treat_g_as_a: bool = False) -> AnalysisResult:
        """Find top N most frequent variants without ambiguities"""
        result = AnalysisResult()
        result.total_sequences = len(self.sequences)
        
        if not self.sequences:
            result.message = "No sequences to analyze"
            return result
        
        counts = Counter(self.sequences)
        total = len(self.sequences)
        
        # Get top N variants
        top_variants = counts.most_common(n)
        
        covered = 0
        for seq, count in top_variants:
            percentage = (count / total) * 100
            result.variants.append((seq, count, percentage))
            covered += count
        
        result.coverage = (covered / total) * 100
        result.uncovered_count = total - covered
        result.uncovered_percentage = 100 - result.coverage
        
        return result


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (supports IUPAC codes)"""
    return ''.join(COMPLEMENT.get(base, base) for base in reversed(seq))


def format_sequence(seq: str, group_size: int = 3) -> str:
    """Format sequence with spacing for readability"""
    return ' '.join(seq[i:i+group_size] for i in range(0, len(seq), group_size))


def format_results_table(result: AnalysisResult, add_spacer: bool = True, 
                         show_percentages: bool = True) -> str:
    """Format analysis results as a table string"""
    lines = []
    
    for i, (seq, count, pct) in enumerate(result.variants, 1):
        formatted_seq = format_sequence(seq) if add_spacer else seq
        if show_percentages:
            lines.append(f"Variant {i}\t{formatted_seq}\t{count:,}\t{pct:.1f}%")
        else:
            lines.append(f"Variant {i}\t{formatted_seq}\t{count:,}")
    
    return '\n'.join(lines)


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================
if __name__ == "__main__":
    # Simple test when run directly
    test_fasta = """
>Seq1
AAAAAAAAA
>Seq2
ATAAAAAAA
>Seq3
ATAAAAAAA
>Seq4
AAAAAAAAA
>Seq5
AAAAAAAAA
>Seq6
AAAAGAAAA
"""
    
    analyzer = OligoAnalyzer()
    report = analyzer.load_sequences(test_fasta)
    
    print("=== Quality Report ===")
    print(f"Original: {report.original_count}")
    print(f"Valid: {report.valid_count}")
    print(f"Majority length: {report.majority_length}")
    print()
    
    print("=== All Variants ===")
    result = analyzer.count_variants()
    print(format_results_table(result))
    print()
    
    print("=== Min Variants (max 1 ambiguity) ===")
    result = analyzer.find_minimum_variants_greedy(1)
    print(format_results_table(result))
    print()
    
    print("=== Min Variants with G-A Wobble ===")
    result = analyzer.find_minimum_variants_greedy(1, treat_g_as_a=True)
    print(format_results_table(result))
    print()
    
    print("=== Top 2 Variants ===")
    result = analyzer.find_top_n_variants(2)
    print(format_results_table(result))
    print(f"Uncovered: {result.uncovered_count} ({result.uncovered_percentage:.1f}%)")
