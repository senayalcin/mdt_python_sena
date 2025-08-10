#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, gzip, difflib, pandas as pd, os, sys

def read_lines(path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rt') as f:
            return f.readlines()
    with open(path, 'r') as f:
        return f.readlines()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--expected', required=True, help='TSV/TSV.GZ/XLSX')
    ap.add_argument('--actual', required=True, help='TSV/TSV.GZ')
    args = ap.parse_args()

    # Load lines (convert xlsx to tsv in-memory if needed)
    if args.expected.endswith('.xlsx'):
        df = pd.read_excel(args.expected)
        exp_tmp = '/mnt/data/_expected_tmp.tsv'
        df.to_csv(exp_tmp, sep='\t', index=False)
        exp_lines = read_lines(exp_tmp)
    else:
        exp_lines = read_lines(args.expected)

    act_lines = read_lines(args.actual)

    print("Expected lines:", len(exp_lines))
    print("Actual lines:", len(act_lines))

    # Save detailed diff
    diff_path = 'output_test/detailed_diff.txt'
    with open(diff_path, 'w') as f:
        diff = difflib.unified_diff(exp_lines, act_lines, fromfile='expected', tofile='actual', lineterm='')
        f.write('\n'.join(diff))
    print(f"Detailed diff saved to {diff_path}")

    # Dataframe summary
    try:
        if args.expected.endswith('.xlsx'):
            exp_df = pd.read_excel(args.expected)
        elif args.expected.endswith('.gz'):
            exp_df = pd.read_csv(args.expected, sep='\t', compression='gzip', comment='#')
        else:
            exp_df = pd.read_csv(args.expected, sep='\t', comment='#')

        if args.actual.endswith('.gz'):
            act_df = pd.read_csv(args.actual, sep='\t', compression='gzip', comment='#')
        else:
            act_df = pd.read_csv(args.actual, sep='\t', comment='#')

        print("Expected shape:", exp_df.shape)
        print("Actual shape:", act_df.shape)
        same_cols = list(exp_df.columns) == list(act_df.columns)
        print("Columns identical:", same_cols)
        if not same_cols:
            print("Expected columns (first 5):", list(exp_df.columns)[:5])
            print("Actual columns (first 5):", list(act_df.columns)[:5])
    except Exception as e:
        print("Dataframe comparison failed:", e, file=sys.stderr)

if __name__ == '__main__':
    main()
