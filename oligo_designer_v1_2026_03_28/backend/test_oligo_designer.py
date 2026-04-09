"""
Tests for oligo_designer.py

The key test: design oligos from a sequence, then reassemble them in silico
and verify you get back the original sequence.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from oligo_designer import (
    tile_sequence, reverse_complement, design_oligos,
    find_homopolymers, find_g_quadruplex, find_repeats,
    gc_content, check_overlap, Overlap, scan_sequence_complexity,
)


def assemble_in_silico(oligos, full_seq):
    """
    Simulate assembly: given a list of Oligo objects, reconstruct the
    full sequence by finding overlaps between adjacent oligos and joining.

    This mimics overlap extension assembly: oligos anneal at overlaps,
    polymerase fills the gaps, producing the full dsDNA product.

    Returns the assembled sequence.
    """
    # Convert all oligos to sense-strand coordinates
    segments = []
    for oligo in oligos:
        if oligo.strand == "antisense":
            seq = reverse_complement(oligo.seq)
        else:
            seq = oligo.seq
        segments.append((oligo.start, oligo.end, seq))

    # Sort by start position
    segments.sort(key=lambda x: x[0])

    # Build consensus by overlaying segments
    assembled = [''] * (max(s[1] for s in segments))
    for start, end, seq in segments:
        for i, base in enumerate(seq):
            pos = start + i
            if pos < len(assembled):
                if assembled[pos] == '':
                    assembled[pos] = base
                else:
                    # Overlap region — should agree
                    assert assembled[pos] == base, (
                        f"Mismatch at position {pos}: "
                        f"existing={assembled[pos]}, new={base}"
                    )

    return ''.join(assembled)


# ---------------------------------------------------------------------------
# Test sequences
# ---------------------------------------------------------------------------

# A real-ish codon-optimized sequence (~500bp)
TEST_SEQ_500 = (
    "ATGCCTGTGCAGAAATCTGTAGCTGAAGCTAAGCCTCCTGAACACTCTGGTGGAGAAGCTGTTTATCCT"
    "GAAGGTGCTGAGTTCGATGCTGCTACCATCGACGAGAAAGCTGATGAAGGTCAGTATATTCAGCCTCTG"
    "CGTGCTTCAGCTGCTCAGGCTCGTATCGCTCGTATTGGTGAAGGTGATGCTATCACTAAAGACGAGACT"
    "GCTGAAATTATCAAGCCTGGTATCAAACGTCTGGATAAACTGCTGGAAGCTAAATCTCCTGAAATCGACA"
    "AAGACGGTATCATCAAAGAGATCAAAGCTCAGATCGATAAAGAAGGTATCCTGACTGCTGAAGAAGTTTC"
    "TAAGATGGCTACCACTGGTAAACGTATCAACGATGCTGGCCGTAGCGCTCTGATCAACGAAGCTATCCGT"
    "CAGATCGAAACCGAAATGGCTATCTCTGCTAACCAGACTCAGGAAGCTGAGCGTAACCTGAAGCAGTCTC"
    "TGATCGATCAGCAGGGTATCGTGAAGACCATCGACATCGACCTG"
)

# Short sequence for edge case testing
TEST_SEQ_SHORT = "ATGCCTGTGCAGAAATCTGTAGCTGAAGCTAAGCCTCCTGAACACTCTGGTGGAGAAGCTGTTTATCCTGAAGGTGCTGAG"

# Sequence with a known problematic region (poly-G)
TEST_SEQ_POLYG = (
    "ATGCCTGTGCAGAAATCTGTAGCTGAAGCTAAGCCTCCTGAACAC"
    "GGGGGGGGG"  # 9x G run - should be flagged
    "TCTGGTGGAGAAGCTGTTTATCCTGAAGGTGCTGAGTTCGATGCT"
    "GCTACCATCGACGAGAAAGCTGATGAAGGTCAGTATATTCAGCCTCTGCGT"
)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_reverse_complement():
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GCTA") == "TAGC"
    print("  PASS: reverse_complement")


def test_gc_content():
    assert gc_content("GGCC") == 1.0
    assert gc_content("AATT") == 0.0
    assert gc_content("ATGC") == 0.5
    assert gc_content("") == 0
    print("  PASS: gc_content")


def test_homopolymers():
    runs = find_homopolymers("AAAAATTTGGGGGGCC", min_run=4)
    assert len(runs) == 2  # AAAAA and GGGGGG
    assert runs[0] == (0, 5, 'A')
    assert runs[1] == (8, 14, 'G')
    print("  PASS: find_homopolymers")


def test_g_quadruplex():
    # Classic G-quadruplex motif
    quads = find_g_quadruplex("GGGATGGGCTGGGATGGG")
    assert len(quads) > 0, "Should detect G-quadruplex"
    # No quadruplex
    quads2 = find_g_quadruplex("ATCGATCGATCG")
    assert len(quads2) == 0
    print("  PASS: find_g_quadruplex")


def test_find_repeats():
    seq = "ATCGATCGATCGATCG"  # "ATCGATCGAT" repeats
    reps = find_repeats(seq, kmer_size=10)
    assert len(reps) > 0, "Should find repeated 10-mers"
    print("  PASS: find_repeats")


def test_in_silico_assembly_500bp():
    """Core test: tile a 500bp sequence, reassemble, verify match."""
    oligos, overlaps = tile_sequence(TEST_SEQ_500)
    assembled = assemble_in_silico(oligos, TEST_SEQ_500)
    assert assembled == TEST_SEQ_500, (
        f"Assembly mismatch!\n"
        f"Expected length: {len(TEST_SEQ_500)}\n"
        f"Got length:      {len(assembled)}\n"
        f"First diff at:   {next(i for i, (a, b) in enumerate(zip(assembled, TEST_SEQ_500)) if a != b) if assembled != TEST_SEQ_500 else 'N/A'}"
    )
    print(f"  PASS: in_silico_assembly_500bp ({len(oligos)} oligos, {len(overlaps)} overlaps)")


def test_in_silico_assembly_short():
    """Short sequence edge case."""
    oligos, overlaps = tile_sequence(TEST_SEQ_SHORT)
    assembled = assemble_in_silico(oligos, TEST_SEQ_SHORT)
    assert assembled == TEST_SEQ_SHORT
    print(f"  PASS: in_silico_assembly_short ({len(oligos)} oligos)")


def test_in_silico_assembly_with_flanks():
    """Test with plasmid flanking sequences."""
    upstream = "GCTAGCTAGCTAGCTAGCTA"  # 20bp
    downstream = "TTAAGCTTGCATGCCTGCAG"  # 20bp
    oligos, overlaps = tile_sequence(
        TEST_SEQ_500,
        plasmid_upstream=upstream,
        plasmid_downstream=downstream,
    )
    full_expected = upstream + TEST_SEQ_500 + downstream
    assembled = assemble_in_silico(oligos, full_expected)
    assert assembled == full_expected, "Assembly with flanks failed!"
    # First oligo should start with upstream sequence
    first_seq = oligos[0].seq  # sense strand
    assert first_seq.startswith(upstream[:10]), "First oligo should include upstream flank"
    print(f"  PASS: in_silico_assembly_with_flanks ({len(oligos)} oligos)")


def test_polyg_detection():
    """Sequence with poly-G should flag issues."""
    issues = scan_sequence_complexity(TEST_SEQ_POLYG)
    homopolymer_issues = [i for i in issues if i["type"] == "homopolymer"]
    assert len(homopolymer_issues) > 0, "Should detect poly-G run"
    print(f"  PASS: polyg_detection ({len(issues)} total issues)")


def test_cross_hybridization():
    """Overlaps with similar sequences should be flagged via thermodynamic check."""
    # Create a repetitive sequence where overlaps will be similar
    repeat = "AGCTGATCGATCGATCGATC"
    spacer = "ATCGATCGATCGATCGATCGATCG"
    seq = repeat + spacer + repeat + spacer + repeat + spacer
    oligos, overlaps = tile_sequence(seq)
    cross_issues = [
        iss for ovl in overlaps for iss in ovl.issues
        if iss.kind == "cross_hybridization"
    ]
    assert len(cross_issues) > 0, "Should detect cross-hybridization in repetitive sequence"
    print(f"  PASS: cross_hybridization ({len(cross_issues)} issues found)")


def test_different_oligo_lengths():
    """Test with 60bp oligos (the longer IDT option)."""
    oligos, overlaps = tile_sequence(TEST_SEQ_500, oligo_length=60, overlap_length=20)
    assembled = assemble_in_silico(oligos, TEST_SEQ_500)
    assert assembled == TEST_SEQ_500
    assert all(o.length <= 80 for o in oligos), "No oligo should be unreasonably long"
    # Should need fewer oligos with longer oligos
    oligos_short, _ = tile_sequence(TEST_SEQ_500, oligo_length=45, overlap_length=20)
    assert len(oligos) < len(oligos_short), "60bp oligos should need fewer total"
    print(f"  PASS: different_oligo_lengths (60bp: {len(oligos)} oligos vs 45bp: {len(oligos_short)} oligos)")


def test_all_overlaps_correct_length():
    """All overlaps should be ~20bp."""
    oligos, overlaps = tile_sequence(TEST_SEQ_500)
    for ovl in overlaps:
        assert 15 <= len(ovl.seq) <= 25, f"Overlap length {len(ovl.seq)} out of range"
    print(f"  PASS: all_overlaps_correct_length")


if __name__ == "__main__":
    print("Running oligo designer tests...\n")
    test_reverse_complement()
    test_gc_content()
    test_homopolymers()
    test_g_quadruplex()
    test_find_repeats()
    test_in_silico_assembly_500bp()
    test_in_silico_assembly_short()
    test_in_silico_assembly_with_flanks()
    test_polyg_detection()
    test_cross_hybridization()
    test_different_oligo_lengths()
    test_all_overlaps_correct_length()
    print("\nAll tests passed!")
