#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script to compare Python output with expected Perl output
"""

import sys
import os
import gzip
import difflib
import pandas as pd
import argparse


def compare_outputs(expected_file, python_output_file):
    print("=" * 80)
    print("OUTPUT COMPARISON ANALYSIS")
    print("=" * 80)

    # 1) Read expected output
    print("\n1. Reading expected output file...")
    try:
        if expected_file.endswith('.gz'):
            with gzip.open(expected_file, 'rt') as f:
                expected_lines = f.readlines()
        else:
            with open(expected_file, 'r') as f:
                expected_lines = f.readlines()
        print(f"   ✓ Expected file has {len(expected_lines)} lines")
    except Exception as e:
        print(f"   ✗ ERROR: Could not read expected file: {expected_file}")
        print(f"     Error: {e}")
        return

    # 2) Read Python output
    print("\n2. Reading Python output file...")
    try:
        if python_output_file.endswith('.gz'):
            with gzip.open(python_output_file, 'rt') as f:
                python_lines = f.readlines()
        else:
            with open(python_output_file, 'r') as f:
                python_lines = f.readlines()
        print(f"   ✓ Python output file has {len(python_lines)} lines")
    except Exception as e:
        print(f"   ✗ ERROR: Could not read Python output file: {python_output_file}")
        print(f"     Error: {e}")
        return

    # 3) Line count comparison
    print("\n3. Line Count Comparison:")
    print(f"   - Expected lines: {len(expected_lines)}")
    print(f"   - Python lines: {len(python_lines)}")
    diff_lines = len(python_lines) - len(expected_lines)
    if diff_lines == 0:
        print(f"   ✓ Line counts match perfectly!")
    else:
        print(f"   ✗ Difference: {diff_lines:+d} lines")

    # 4) Header comparison
    print("\n4. Header Comparison:")
    if expected_lines and python_lines:
        if expected_lines[0] == python_lines[0]:
            print("   ✓ Headers match perfectly")
        else:
            print("   ✗ Headers differ")
            print(f"   Expected header (first 150 chars):\n     {expected_lines[0][:150]}...")
            print(f"   Python header (first 150 chars):\n     {python_lines[0][:150]}...")

    # 5) Content analysis (excluding header)
    print("\n5. Content Analysis:")
    if len(expected_lines) > 1 and len(python_lines) > 1:
        expected_content = expected_lines[1:]
        python_content = python_lines[1:]

        exact_matches = 0
        differences = []

        # Compare line by line
        for i in range(min(len(expected_content), len(python_content))):
            if expected_content[i] == python_content[i]:
                exact_matches += 1
            else:
                differences.append((i + 2, expected_content[i], python_content[i]))

        total_compared = min(len(expected_content), len(python_content))
        match_percentage = (exact_matches / total_compared * 100) if total_compared > 0 else 0

        print(f"   - Exact line matches: {exact_matches}/{total_compared} ({match_percentage:.2f}%)")
        print(f"   - Lines with differences: {len(differences)}")

        if differences:
            print("\n6. First 10 Differences Found:")
            for i, (line_num, exp_line, py_line) in enumerate(differences[:10]):
                print(f"\n   Difference #{i + 1} at line {line_num}:")
                exp_fields = exp_line.strip().split('\t')
                py_fields = py_line.strip().split('\t')
                print(f"   - Expected has {len(exp_fields)} fields")
                print(f"   - Python has {len(py_fields)} fields")

                # Find which fields differ
                diff_fields = []
                for j in range(min(len(exp_fields), len(py_fields))):
                    if exp_fields[j] != py_fields[j]:
                        diff_fields.append(j + 1)

                if diff_fields:
                    print(f"   - Fields that differ: {diff_fields[:5]}")
                    for j in diff_fields[:3]:  # Show first 3 field differences
                        if j - 1 < len(exp_fields) and j - 1 < len(py_fields):
                            print(f"     Field {j}:")
                            print(f"       Expected: {exp_fields[j - 1][:60]}")
                            print(f"       Python:   {py_fields[j - 1][:60]}")

    # 7) Statistical summary
    print("\n7. Statistical Summary:")
    try:
        # Read as dataframes
        if expected_file.endswith('.gz'):
            expected_df = pd.read_csv(expected_file, sep='\t', comment='#', compression='gzip')
        else:
            expected_df = pd.read_csv(expected_file, sep='\t', comment='#')

        if python_output_file.endswith('.gz'):
            python_df = pd.read_csv(python_output_file, sep='\t', comment='#', compression='gzip')
        else:
            python_df = pd.read_csv(python_output_file, sep='\t', comment='#')

        print(f"   - Expected data shape: {expected_df.shape} (rows × columns)")
        print(f"   - Python data shape: {python_df.shape} (rows × columns)")

        # Column comparison
        if list(expected_df.columns) == list(python_df.columns):
            print("   ✓ Column names match perfectly")
        else:
            print("   ✗ Column names differ")
            exp_cols = set(expected_df.columns)
            py_cols = set(python_df.columns)
            missing_cols = exp_cols - py_cols
            extra_cols = py_cols - exp_cols
            if missing_cols:
                print(f"     Missing columns: {list(missing_cols)[:5]}")
            if extra_cols:
                print(f"     Extra columns: {list(extra_cols)[:5]}")

        # First column (key) analysis
        if expected_df.shape[1] > 0 and python_df.shape[1] > 0:
            key_col = expected_df.columns[0]
            if key_col in python_df.columns:
                expected_unique = set(expected_df[key_col].unique())
                python_unique = set(python_df[key_col].unique())
                common = expected_unique & python_unique
                only_expected = expected_unique - python_unique
                only_python = python_unique - expected_unique

                print(f"\n   Key column '{key_col}' analysis:")
                print(f"   - Common values: {len(common)}")
                print(f"   - Only in expected: {len(only_expected)}")
                print(f"   - Only in Python: {len(only_python)}")

                if only_expected and len(only_expected) <= 10:
                    print(f"   - Missing from Python: {sorted(list(only_expected))[:10]}")
                elif only_expected:
                    print(f"   - First 10 missing from Python: {sorted(list(only_expected))[:10]}")

                if only_python and len(only_python) <= 10:
                    print(f"   - Extra in Python: {sorted(list(only_python))[:10]}")
                elif only_python:
                    print(f"   - First 10 extra in Python: {sorted(list(only_python))[:10]}")

    except Exception as e:
        print(f"   ✗ Could not perform statistical analysis: {e}")

    # 8) Generate detailed diff file
    print("\n8. Generating Detailed Diff File...")

    # Create output directory if it doesn't exist
    os.makedirs('output_test', exist_ok=True)

    diff_file = 'output_test/detailed_diff.txt'
    try:
        with open(diff_file, 'w') as f:
            diff = difflib.unified_diff(
                expected_lines,
                python_lines,
                fromfile='expected_output',
                tofile='python_output',
                lineterm='',
                n=3  # Context lines
            )
            diff_content = '\n'.join(diff)
            f.write(diff_content)

        # Check file size
        file_size = os.path.getsize(diff_file)
        if file_size > 0:
            print(f"   ✓ Detailed diff saved to '{diff_file}' ({file_size:,} bytes)")
        else:
            print(f"   ⚠ Diff file created but is empty (files might be identical)")

    except Exception as e:
        print(f"   ✗ Could not create diff file: {e}")

    # 9) Summary
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)

    if len(expected_lines) == len(python_lines) and exact_matches == total_compared:
        print("✓ PERFECT MATCH! Output files are identical.")
    elif match_percentage >= 99:
        print(f"✓ EXCELLENT! {match_percentage:.2f}% of lines match.")
    elif match_percentage >= 95:
        print(f"⚠ GOOD: {match_percentage:.2f}% of lines match.")
    else:
        print(f"✗ NEEDS WORK: Only {match_percentage:.2f}% of lines match.")


def main():
    parser = argparse.ArgumentParser(
        description="Compare Python output with expected Perl output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python test_comparison.py --expected perl_output.tsv.gz --python python_output.tsv.gz
  python test_comparison.py --expected expected.tsv --python actual.tsv
        """
    )

    parser.add_argument('--expected', required=True,
                        help='Path to expected output file (TSV or TSV.GZ)')
    parser.add_argument('--python', required=True,
                        help='Path to Python output file (TSV or TSV.GZ)')
    parser.add_argument('--verbose', action='store_true',
                        help='Show more detailed output')

    args = parser.parse_args()

    # Check if files exist
    if not os.path.exists(args.expected):
        print(f"Error: Expected file '{args.expected}' does not exist!")
        sys.exit(1)

    if not os.path.exists(args.python):
        print(f"Error: Python output file '{args.python}' does not exist!")
        sys.exit(1)

    # Run comparison
    compare_outputs(args.expected, args.python)

    # Print example commands
    print("\n" + "=" * 80)
    print("NEXT STEPS")
    print("=" * 80)
    print("\nTo generate Python output with the same parameters as Perl:")
    print("""
python mostDominantTranscriptSwitch.py \\
  --pcawg Datas/pcawg.rnaseq.transcript.expr.tpm.tsv.gz \\
  --gtex Datas/GTEX_v4.pcawg.transcripts.tpm.tsv.gz \\
  --ensg Datas/ensg_ensp_enst_ense_geneName_v75.tsv.gz \\
  --canon Datas/canonEnsp_ensg_ensp_enst_geneName_v75.tsv.gz \\
  --iso Datas/interactionsInIsoforms_900_2.tsv.gz \\
  --seq Datas/ensp_ensg_enst_sequence.tsv.gz \\
  --redundant Datas/redundantENST_v75.txt \\
  --minString 900 \\
  --minEnr 2 \\
  --maxQ 0.01 \\
  --minExp 2 \\
  --enst 1 \\
  | gzip > output_py/python_output.tsv.gz
""")

    print("\nThen compare outputs:")
    print("python test_comparison.py --expected perl_output.tsv.gz --python output_py/python_output.tsv.gz")


if __name__ == "__main__":
    main()