#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script to compare Python output with expected output
"""

import sys
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
        print(f"   - Expected file has {len(expected_lines)} lines")
    except Exception as e:
        print(f"   ERROR: Could not read expected file: {expected_file} ({e})")
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
        print(f"   - Python output file has {len(python_lines)} lines")
    except Exception as e:
        print(f"   ERROR: Could not read Python output file: {python_output_file} ({e})")
        return

    # 3) Line count comparison
    print("\n3. Line Count Comparison:")
    print(f"   - Expected lines: {len(expected_lines)}")
    print(f"   - Python lines: {len(python_lines)}")
    print(f"   - Difference: {len(python_lines) - len(expected_lines)} lines")

    # 4) Header comparison
    print("\n4. Header Comparison:")
    if expected_lines and python_lines:
        if expected_lines[0] == python_lines[0]:
            print("   ✓ Headers match perfectly")
        else:
            print("   ✗ Headers differ")
            print(f"   Expected header: {expected_lines[0][:100]}...")
            print(f"   Python header: {python_lines[0][:100]}...")

    # 5) Content analysis (excluding header)
    print("\n5. Content Analysis:")
    if len(expected_lines) > 1 and len(python_lines) > 1:
        expected_content = expected_lines[1:]
        python_content = python_lines[1:]

        exact_matches = 0
        differences = []
        for i in range(min(len(expected_content), len(python_content))):
            if expected_content[i] == python_content[i]:
                exact_matches += 1
            else:
                differences.append((i + 2, expected_content[i], python_content[i]))

        print(f"   - Exact line matches: {exact_matches}/{min(len(expected_content), len(python_content))}")
        print(f"   - Lines with differences: {len(differences)}")

        if differences:
            print("\n6. First 5 Differences Found:")
            for i, (line_num, exp_line, py_line) in enumerate(differences[:5]):
                print(f"\n   Difference #{i + 1} at line {line_num}:")
                exp_fields = exp_line.strip().split('\t')
                py_fields = py_line.strip().split('\t')
                print(f"   Expected fields: {len(exp_fields)}")
                print(f"   Python fields: {len(py_fields)}")
                for j in range(min(len(exp_fields), len(py_fields))):
                    if exp_fields[j] != py_fields[j]:
                        print(f"   Field {j + 1} differs:")
                        print(f"     Expected: {exp_fields[j][:50]}...")
                        print(f"     Python:   {py_fields[j][:50]}...")

    # 7) Statistical summary
    print("\n7. Statistical Summary:")
    try:
        expected_df = pd.read_csv(expected_file, sep='\t', comment='#')  # pandas infers .gz
        if python_output_file.endswith('.gz'):
            python_df = pd.read_csv(python_output_file, sep='\t', comment='#', compression='gzip')
        else:
            python_df = pd.read_csv(python_output_file, sep='\t', comment='#')

        print(f"   - Expected data shape: {expected_df.shape}")
        print(f"   - Python data shape: {python_df.shape}")

        if list(expected_df.columns) == list(python_df.columns):
            print("   ✓ Column names match")
        else:
            print("   ✗ Column names differ")
            print(f"     Expected columns: {list(expected_df.columns)[:5]}...")
            print(f"     Python columns: {list(python_df.columns)[:5]}...")

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
                if only_expected:
                    print(f"   - First 5 missing from Python: {list(only_expected)[:5]}")
                if only_python:
                    print(f"   - First 5 extra in Python: {list(only_python)[:5]}")
    except Exception as e:
        print(f"   Could not perform statistical analysis: {e}")

    # 8) Detailed diff
    print("\n8. Generating Detailed Diff File...")
    with open('output_test/detailed_diff.txt', 'w') as f:
        diff = difflib.unified_diff(expected_lines, python_lines,
                                    fromfile='expected_output',
                                    tofile='python_output',
                                    lineterm='')
        f.write('\n'.join(diff))
        print("   - Detailed diff saved to 'detailed_diff.txt'")

    print("\n" + "=" * 80)
    print("COMPARISON COMPLETE")
    print("=" * 80)

def main():
    ap = argparse.ArgumentParser(description="Compare expected output with main11.py output")
    ap.add_argument("--expected", required=True, help="Path to expected TSV/TSV.GZ")
    ap.add_argument("--python-output", required=True, help="Path to main11.py output (TSV or TSV.GZ)")
    args = ap.parse_args()

    compare_outputs(args.expected, args.python_output)

    print("\n\nValidation Script Commands:")
    print("-" * 80)
    print("To run the Python script with the same parameters as Perl:")
    print(r"""
python main11.py \
  --pcawg pcawg.rnaseq.transcript.expr.tpm.tsv.gz \
  --gtex GTEX_v4.pcawg.transcripts.tpm.tsv.gz \
  --ensg ensg_ensp_enst_ense_geneName_v75.tsv.gz \
  --canon canonEnsp_ensg_ensp_enst_geneName_v75.tsv.gz \
  --iso interactionsInIsoforms_900_2.tsv.gz \
  --seq ensp_ensg_enst_sequence.tsv.gz \
  --redundant redundantENST_v75.txt \
  --minString 900 \
  --minEnr 2 \
  --maxQ 0.01 \
  --minExp 2 \
  --enst 1 \
  --verbose \
  | gzip > output_test/python_output.tsv.gz
""")
    print("Then run this comparison script:")
    print("python main11_test.py --expected <EXPECTED.tsv[.gz]> --python-output python_output.tsv.gz")

if __name__ == "__main__":
    main()
