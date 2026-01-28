#!/usr/bin/env python3
"""
Comprehensive tests for the Oligo Analyzer Engine
"""

import random
from oligo_engine import OligoAnalyzer, format_sequence, format_results_table


def generate_random_sequences(base_seq: str, n_sequences: int, mutation_rate: float = 0.1) -> str:
    """Generate random test sequences based on a base sequence"""
    bases = ['A', 'C', 'G', 'T']
    fasta_lines = []
    
    for i in range(n_sequences):
        seq = list(base_seq)
        for j in range(len(seq)):
            if random.random() < mutation_rate:
                seq[j] = random.choice(bases)
        fasta_lines.append(f">Seq{i+1}")
        fasta_lines.append(''.join(seq))
    
    return '\n'.join(fasta_lines)


def test_basic_functionality():
    """Test basic counting functionality"""
    print("=" * 60)
    print("TEST 1: Basic Functionality")
    print("=" * 60)
    
    test_data = """
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
    report = analyzer.load_sequences(test_data)
    
    print(f"Loaded {report.valid_count} valid sequences")
    assert report.valid_count == 6, "Expected 6 valid sequences"
    
    # Test all variants
    result = analyzer.count_variants()
    print(f"Found {len(result.variants)} unique variants")
    assert len(result.variants) == 3, "Expected 3 unique variants"
    
    print("✓ Basic functionality test PASSED\n")


def test_quality_filtering():
    """Test quality filtering with various issues"""
    print("=" * 60)
    print("TEST 2: Quality Filtering")
    print("=" * 60)
    
    test_data = """
>Good1
AAAAAAAAA
>Good2
AAAAAAAAA
>WithGap
AAA-AAAAA
>WithAmbiguity
AAANAAAAA
>WrongLength
AAAA
>AnotherGood
AAAAAAAAA
"""
    
    analyzer = OligoAnalyzer()
    report = analyzer.load_sequences(test_data)
    
    print(f"Original: {report.original_count}")
    print(f"Valid: {report.valid_count}")
    print(f"Removed gaps: {report.removed_gaps}")
    print(f"Removed ambiguous: {report.removed_ambiguous}")
    print(f"Removed wrong length: {report.removed_wrong_length}")
    
    assert report.original_count == 6, "Expected 6 original sequences"
    assert report.valid_count == 3, "Expected 3 valid sequences"
    assert report.removed_gaps == 1, "Expected 1 sequence removed for gaps"
    assert report.removed_ambiguous == 1, "Expected 1 sequence removed for ambiguity"
    assert report.removed_wrong_length == 1, "Expected 1 sequence removed for wrong length"
    
    print("✓ Quality filtering test PASSED\n")


def test_ambiguity_consolidation():
    """Test ambiguity-based consolidation"""
    print("=" * 60)
    print("TEST 3: Ambiguity Consolidation")
    print("=" * 60)
    
    test_data = """
>Seq1
AAAAAAAAA
>Seq2
ATAAAAAAA
>Seq3
ACAAAAAAA
>Seq4
AGAAAAAAA
"""
    
    analyzer = OligoAnalyzer()
    analyzer.load_sequences(test_data)
    
    # With 1 ambiguity, should consolidate all 4 into one with N at position 2
    result = analyzer.find_minimum_variants_greedy(1)
    print(f"With max 1 ambiguity: {len(result.variants)} variant(s)")
    print(format_results_table(result))
    
    # With 0 ambiguities, should have 4 variants
    result = analyzer.find_minimum_variants_greedy(0)
    print(f"\nWith max 0 ambiguities: {len(result.variants)} variants")
    assert len(result.variants) == 4, "Expected 4 variants with 0 ambiguities"
    
    print("✓ Ambiguity consolidation test PASSED\n")


def test_g_a_wobble():
    """Test G-A wobble base pairing option"""
    print("=" * 60)
    print("TEST 4: G-A Wobble Base Pairing")
    print("=" * 60)
    
    test_data = """
>Seq1
AAAAAAAAA
>Seq2
GAAAAAAAA
>Seq3
AAAGAAAAA
>Seq4
AAAAAAGAA
"""
    
    analyzer = OligoAnalyzer()
    analyzer.load_sequences(test_data)
    
    # Without G-A wobble
    result = analyzer.find_minimum_variants_greedy(0, treat_g_as_a=False)
    print(f"Without G-A wobble: {len(result.variants)} variants")
    
    # With G-A wobble
    result = analyzer.find_minimum_variants_greedy(0, treat_g_as_a=True)
    print(f"With G-A wobble: {len(result.variants)} variant(s)")
    assert len(result.variants) == 1, "Expected 1 variant with G-A wobble"
    
    print("✓ G-A wobble test PASSED\n")


def test_top_n():
    """Test top N variants functionality"""
    print("=" * 60)
    print("TEST 5: Top N Variants")
    print("=" * 60)
    
    test_data = """
>Seq1
AAAAAAAAA
>Seq2
AAAAAAAAA
>Seq3
AAAAAAAAA
>Seq4
TTTTTTTTT
>Seq5
TTTTTTTTT
>Seq6
CCCCCCCCC
"""
    
    analyzer = OligoAnalyzer()
    analyzer.load_sequences(test_data)
    
    result = analyzer.find_top_n_variants(2)
    print(f"Top 2 variants cover {result.coverage:.1f}% of sequences")
    print(f"Uncovered: {result.uncovered_count} ({result.uncovered_percentage:.1f}%)")
    print(format_results_table(result))
    
    assert len(result.variants) == 2, "Expected 2 variants"
    assert result.uncovered_count == 1, "Expected 1 uncovered sequence"
    
    print("✓ Top N test PASSED\n")


def test_large_dataset():
    """Test performance with larger dataset"""
    print("=" * 60)
    print("TEST 6: Large Dataset Performance")
    print("=" * 60)
    
    import time
    
    # Generate 500 sequences with some variation
    random.seed(42)  # For reproducibility
    base_seq = "ACGTACGTACGTACGTACGT"  # 20bp
    test_data = generate_random_sequences(base_seq, 500, mutation_rate=0.05)
    
    analyzer = OligoAnalyzer()
    
    start = time.time()
    report = analyzer.load_sequences(test_data)
    load_time = time.time() - start
    print(f"Loaded {report.valid_count} sequences in {load_time:.3f}s")
    
    start = time.time()
    result = analyzer.count_variants()
    count_time = time.time() - start
    print(f"Counted {len(result.variants)} unique variants in {count_time:.3f}s")
    
    start = time.time()
    result = analyzer.find_minimum_variants_greedy(2)
    greedy_time = time.time() - start
    print(f"Consolidated to {len(result.variants)} variants (max 2 amb) in {greedy_time:.3f}s")
    
    print("✓ Large dataset test PASSED\n")


def test_format_functions():
    """Test formatting functions"""
    print("=" * 60)
    print("TEST 7: Formatting Functions")
    print("=" * 60)
    
    seq = "ACGTACGTACGT"
    formatted = format_sequence(seq, 3)
    print(f"Original: {seq}")
    print(f"Formatted: {formatted}")
    assert formatted == "ACG TAC GTA CGT", "Formatting incorrect"
    
    print("✓ Formatting test PASSED\n")


def test_edge_cases():
    """Test edge cases"""
    print("=" * 60)
    print("TEST 8: Edge Cases")
    print("=" * 60)
    
    analyzer = OligoAnalyzer()
    
    # Empty input
    report = analyzer.load_sequences("")
    print(f"Empty input: {report.valid_count} sequences")
    assert report.valid_count == 0
    
    # Single sequence
    report = analyzer.load_sequences(">Seq1\nAAAA")
    print(f"Single sequence: {report.valid_count} sequences")
    assert report.valid_count == 1
    
    result = analyzer.count_variants()
    assert len(result.variants) == 1
    
    # All identical
    report = analyzer.load_sequences(">S1\nAAAA\n>S2\nAAAA\n>S3\nAAAA")
    result = analyzer.count_variants()
    print(f"All identical: {len(result.variants)} variant(s)")
    assert len(result.variants) == 1
    
    # All different
    report = analyzer.load_sequences(">S1\nAAAA\n>S2\nTTTT\n>S3\nCCCC\n>S4\nGGGG")
    result = analyzer.count_variants()
    print(f"All different: {len(result.variants)} variants")
    assert len(result.variants) == 4
    
    print("✓ Edge cases test PASSED\n")


def test_example_from_spec():
    """Test the exact examples from the specification"""
    print("=" * 60)
    print("TEST 9: Specification Examples")
    print("=" * 60)
    
    # Note: Spaces in sequences are ignored
    test_data = """
>Sequence1
AAA AAA AAA
>Sequence2
ATA AAA AAA
>Sequence3
ATA AAA AAA
>Sequence4
AAA AAA AAA
>Sequence5
AAA AAA AAA
>Sequence6
AAA AGA AAA
"""
    
    analyzer = OligoAnalyzer()
    analyzer.load_sequences(test_data)
    
    # Example 1.1: All variations without ambiguities
    print("Example 1.1: All variations without ambiguities")
    result = analyzer.count_variants()
    print(format_results_table(result))
    assert len(result.variants) == 3
    print()
    
    # Example 1.2: Minimum variations with up to 1 ambiguity
    print("Example 1.2: Minimum variations with up to 1 ambiguity")
    result = analyzer.find_minimum_variants_greedy(1)
    print(format_results_table(result))
    assert len(result.variants) == 2
    print()
    
    # Example 1.3: With G-A wobble
    print("Example 1.3: Minimum variations with G-A wobble")
    result = analyzer.find_minimum_variants_greedy(1, treat_g_as_a=True)
    print(format_results_table(result))
    # With G-A wobble and 1 ambiguity, should consolidate to 1 variant
    assert len(result.variants) == 1
    print()
    
    # Example 1.4: Maximum 2 variations
    print("Example 1.4: Top 2 variations")
    result = analyzer.find_top_n_variants(2)
    print(format_results_table(result))
    print(f"Sequences with mismatches: {result.uncovered_count} ({result.uncovered_percentage:.1f}%)")
    assert len(result.variants) == 2
    
    print("✓ Specification examples test PASSED\n")


def run_all_tests():
    """Run all tests"""
    print("\n" + "=" * 60)
    print("OLIGO ANALYZER TEST SUITE")
    print("=" * 60 + "\n")
    
    test_basic_functionality()
    test_quality_filtering()
    test_ambiguity_consolidation()
    test_g_a_wobble()
    test_top_n()
    test_large_dataset()
    test_format_functions()
    test_edge_cases()
    test_example_from_spec()
    
    print("=" * 60)
    print("ALL TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
